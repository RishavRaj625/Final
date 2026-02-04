from flask import Flask, request, jsonify, render_template
from flask_cors import CORS
import pandas as pd
import joblib
import json
import shap
import numpy as np
from Bio.Data import CodonTable
from collections import defaultdict

# ================= APP SETUP =================

app = Flask(__name__, template_folder=".")
CORS(app)

# ================= LOAD CODON USAGE DATA =================

codon_df = pd.read_csv("codon_usage.csv", low_memory=False)
codon_df.columns = [c.strip().upper() for c in codon_df.columns]

# ================= FIX DATA TYPES =================
# Convert all codon columns to numeric, handling any non-numeric values

def fix_data_types(df):
    """Convert all codon columns to numeric types"""
    # Identify non-metadata columns
    metadata_cols = ["SPECIESNAME", "KINGDOM", "SPECIES", "ORGANISM"]
    
    for col in df.columns:
        if col not in metadata_cols:
            # Convert to numeric, replacing any non-numeric values with 0
            df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0)
    
    return df

codon_df = fix_data_types(codon_df)

# ================= LOAD TRAINED MODEL =================

model = joblib.load("model_outputs/global_codon_bwt_model.pkl")

# Required for ml_weight (frontend)
booster = model.get_booster()
feature_importance = booster.get_score(importance_type="weight")

# ================= LOAD EVALUATION METRICS =================

with open("model_outputs/evaluation_metrics.json", "r") as f:
    EVAL_METRICS = json.load(f)

# ================= GENETIC CODE =================

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

table = CodonTable.unambiguous_rna_by_id[1]
AA_TO_CODONS = {}
for codon, aa in table.forward_table.items():
    AA_TO_CODONS.setdefault(aa, []).append(codon)

# ================= FEATURE COLUMNS (FOR SHAP) =================

# More robust feature column detection
metadata_cols = ["SPECIESNAME", "KINGDOM", "SPECIES", "ORGANISM"]
FEATURE_COLUMNS = [
    c for c in codon_df.columns
    if c not in metadata_cols and codon_df[c].dtype in ['float64', 'int64', 'float32', 'int32']
]

print(f"✓ Found {len(FEATURE_COLUMNS)} feature columns")
print(f"✓ Data types verified: {codon_df[FEATURE_COLUMNS].dtypes.unique()}")

# ================= ENHANCED SHAP SETUP =================

# Create stratified background dataset (sample from different species categories)
def create_stratified_background(df, n_samples=100):
    """Create a stratified background dataset for better SHAP explanations"""
    backgrounds = []
    
    # If KINGDOM column exists, stratify by kingdom
    if 'KINGDOM' in df.columns:
        kingdoms = df['KINGDOM'].dropna().unique()
        if len(kingdoms) > 0:
            samples_per_kingdom = max(1, n_samples // len(kingdoms))
            
            for kingdom in kingdoms:
                kingdom_df = df[df['KINGDOM'] == kingdom]
                n = min(samples_per_kingdom, len(kingdom_df))
                if n > 0:
                    backgrounds.append(kingdom_df.sample(n, random_state=42))
            
            if backgrounds:
                X_background = pd.concat(backgrounds)[FEATURE_COLUMNS]
            else:
                X_background = df[FEATURE_COLUMNS].sample(min(n_samples, len(df)), random_state=42)
        else:
            X_background = df[FEATURE_COLUMNS].sample(min(n_samples, len(df)), random_state=42)
    else:
        # Random sampling if no kingdom information
        X_background = df[FEATURE_COLUMNS].sample(min(n_samples, len(df)), random_state=42)
    
    # Ensure all columns are numeric
    X_background = X_background.apply(pd.to_numeric, errors='coerce').fillna(0)
    
    return X_background

print("Creating background dataset...")
X_background = create_stratified_background(codon_df, n_samples=100)
print(f"✓ Background dataset shape: {X_background.shape}")
print(f"✓ Background data types: {X_background.dtypes.unique()}")

print("Initializing SHAP explainer...")
explainer = shap.TreeExplainer(model)
print("✓ SHAP explainer initialized")

# Pre-compute SHAP values for background (for faster queries)
print("Computing background SHAP values...")
try:
    background_shap_values = explainer.shap_values(X_background)
    background_expected_value = explainer.expected_value
    print("✓ Background SHAP values computed")
except Exception as e:
    print(f"⚠ Warning: Could not pre-compute SHAP values: {e}")
    background_shap_values = None
    background_expected_value = None

# Cache for SHAP computations
shap_cache = {}

# ================= CORE CODON ANALYSIS =================

def analyze(aa, selected_codon=None):
    aa = aa.upper()
    selected_codon = selected_codon.upper() if selected_codon else None

    if aa not in VALID_AA:
        return None, "Invalid amino acid (use single-letter code)"

    codons = [c for c in AA_TO_CODONS[aa] if c in codon_df.columns]
    if not codons:
        return None, "No codon data available"

    mean_usage = codon_df[codons].mean().sort_values(ascending=False)

    codon_ranking = []
    selected_rank = None

    for i, codon in enumerate(mean_usage.index, start=1):
        if codon == selected_codon:
            selected_rank = i

        codon_ranking.append({
            "rank": i,
            "codon": codon,
            "frequency": float(mean_usage[codon]),
            "ml_weight": float(
                feature_importance.get(codon.replace("U", "T"), 0.0)
            )
        })

    df_tmp = codon_df.copy()
    df_tmp["SCORE"] = (
        df_tmp[selected_codon]
        if selected_codon in codons
        else df_tmp[codons].sum(axis=1)
    )

    top_species = (
        df_tmp.sort_values("SCORE", ascending=False)
        .head(5)[["SPECIESNAME", "SCORE"]]
        .to_dict(orient="records")
    )

    return {
        "codon_ranking": codon_ranking,
        "selected_rank": selected_rank,
        "top_species": top_species
    }, None

# ================= SPECIES-SPECIFIC PREFERENCE =================

def species_specific_analysis(aa, codon=None):
    codons = [c for c in AA_TO_CODONS.get(aa, []) if c in codon_df.columns]
    if not codons:
        return None

    df = codon_df.copy()
    df["PREFERENCE_SCORE"] = (
        df[codon] if codon in codons else df[codons].sum(axis=1)
    )

    max_val = df["PREFERENCE_SCORE"].max()
    if max_val > 0:
        df["PREFERENCE_SCORE"] /= max_val

    return {
        "top_species": df.sort_values("PREFERENCE_SCORE", ascending=False)
            .head(5)[["SPECIESNAME", "PREFERENCE_SCORE"]]
            .to_dict(orient="records"),
        "bottom_species": df.sort_values("PREFERENCE_SCORE")
            .head(5)[["SPECIESNAME", "PREFERENCE_SCORE"]]
            .to_dict(orient="records"),
        "explanation": "Normalized species-level codon preference analysis."
    }

# ================= HOST-AWARE OPTIMIZATION =================

def host_aware_optimization(aa, host):
    if not host:
        return None

    codons = [c for c in AA_TO_CODONS.get(aa, []) if c in codon_df.columns]
    df = codon_df[codon_df["SPECIESNAME"].str.contains(host, case=False, na=False)]

    if df.empty:
        return None

    mean_usage = df[codons].mean().sort_values(ascending=False)

    return {
        "host_species": host,
        "optimal_codon": mean_usage.index[0],
        "codon_ranking": list(mean_usage.items())
    }

# ================= RARE CODON RISK =================

def rare_codon_risk(aa, host):
    if not host:
        return None

    codons = [c for c in AA_TO_CODONS.get(aa, []) if c in codon_df.columns]
    df = codon_df[codon_df["SPECIESNAME"].str.contains(host, case=False, na=False)]

    if df.empty:
        return None

    mean_usage = df[codons].mean()
    threshold = mean_usage.quantile(0.25)

    rare_codons = mean_usage[mean_usage <= threshold].sort_values()

    return {
        "threshold": float(threshold),
        "rare_codons": [
            {"codon": c, "usage": float(v)}
            for c, v in rare_codons.items()
        ]
    }

# ================= HOST COMPATIBILITY SCORE =================

def host_compatibility_score(aa, host):
    if not host:
        return None

    codons = [c for c in AA_TO_CODONS.get(aa, []) if c in codon_df.columns]
    host_df = codon_df[codon_df["SPECIESNAME"].str.contains(host, case=False, na=False)]

    if host_df.empty:
        return None

    host_profile = host_df[codons].mean()
    global_profile = codon_df[codons].mean()

    similarity = 1 - (abs(host_profile - global_profile).sum() / global_profile.sum())

    return {
        "host_species": host,
        "compatibility_score": round(float(similarity * 100), 2)
    }

# ================= ENHANCED SHAP EXPLANATIONS =================

def shap_explain_codon(codon, species=None, detailed=True):
    """
    Enhanced SHAP explanation for a specific codon
    
    Args:
        codon: RNA codon to explain (e.g., "UUA")
        species: Optional species name for context-specific explanation
        detailed: If True, return detailed explanation with interactions
    
    Returns:
        Dictionary with SHAP values and explanations
    """
    if not codon:
        return None

    feature = codon.replace("U", "T")
    if feature not in FEATURE_COLUMNS:
        return None
    
    # Use cached SHAP values if available
    cache_key = f"{codon}_{species}_{detailed}"
    if cache_key in shap_cache:
        return shap_cache[cache_key]

    try:
        # Select appropriate sample for explanation
        if species:
            # Get species-specific sample
            species_df = codon_df[codon_df["SPECIESNAME"].str.contains(species, case=False, na=False)]
            if not species_df.empty:
                X_sample = species_df[FEATURE_COLUMNS].mean().to_frame().T
            else:
                X_sample = X_background.mean().to_frame().T
        else:
            # Use global average
            X_sample = X_background.mean().to_frame().T
        
        # Ensure sample is numeric
        X_sample = X_sample.apply(pd.to_numeric, errors='coerce').fillna(0)
        
        # Compute SHAP values
        shap_values = explainer.shap_values(X_sample)
        
        # Get contributions for all features
        contributions = dict(zip(FEATURE_COLUMNS, shap_values[0]))
        
        # Sort by absolute contribution
        sorted_features = sorted(
            contributions.items(),
            key=lambda x: abs(x[1]),
            reverse=True
        )
        
        # Use pre-computed expected value if available, otherwise compute
        if background_expected_value is not None:
            expected_val = background_expected_value
        else:
            expected_val = explainer.expected_value
        
        # Base explanation (top 10 features)
        explanation = {
            "codon": codon,
            "base_value": float(expected_val),
            "prediction": float(expected_val + sum(shap_values[0])),
            "total_effect": float(sum(shap_values[0])),
            "top_features": [
                {
                    "feature": f.replace("T", "U"),  # Convert back to RNA
                    "shap_value": round(float(v), 6),
                    "abs_shap_value": round(abs(float(v)), 6),
                    "effect": "positive" if v > 0 else "negative",
                    "percentage": round(abs(float(v)) / (abs(sum(shap_values[0])) + 1e-10) * 100, 2)
                }
                for f, v in sorted_features[:10]
            ]
        }
        
        if detailed:
            # Add interpretation
            explanation["interpretation"] = generate_interpretation(
                codon, sorted_features[:5], species
            )
            
            # Interactions are computationally expensive, so we'll skip them for now
            # or compute them on-demand
            explanation["interactions"] = []
        
        # Cache the result
        shap_cache[cache_key] = explanation
        
        return explanation
    
    except Exception as e:
        print(f"Error computing SHAP explanation: {e}")
        return None


def generate_interpretation(codon, top_features, species=None):
    """
    Generate human-readable interpretation of SHAP values
    """
    if not top_features:
        return "No significant features found."
    
    interpretation = []
    
    # Analyze the top feature
    top_feature, top_value = top_features[0]
    direction = "increases" if top_value > 0 else "decreases"
    
    context = f"for {species}" if species else "globally"
    
    interpretation.append(
        f"The codon {codon} {direction} usage prediction {context}. "
        f"The strongest influence comes from {top_feature.replace('T', 'U')} "
        f"(SHAP value: {top_value:.4f})."
    )
    
    # Analyze positive vs negative contributors
    positive_features = [f for f, v in top_features if v > 0]
    negative_features = [f for f, v in top_features if v < 0]
    
    if positive_features:
        interpretation.append(
            f"Positive contributors: {', '.join([f.replace('T', 'U') for f in positive_features[:3]])}."
        )
    
    if negative_features:
        interpretation.append(
            f"Negative contributors: {', '.join([f.replace('T', 'U') for f in negative_features[:3]])}."
        )
    
    return " ".join(interpretation)


def shap_compare_codons(codons_list, species=None):
    """
    Compare SHAP explanations for multiple codons
    """
    if not codons_list:
        return None
    
    comparisons = []
    
    for codon in codons_list:
        explanation = shap_explain_codon(codon, species, detailed=False)
        if explanation:
            comparisons.append({
                "codon": codon,
                "prediction": explanation["prediction"],
                "total_effect": explanation["total_effect"],
                "top_3_features": explanation["top_features"][:3]
            })
    
    return {
        "comparisons": comparisons,
        "species": species if species else "global",
        "most_influential_codon": max(
            comparisons, 
            key=lambda x: abs(x["total_effect"])
        )["codon"] if comparisons else None
    }


def shap_global_importance(top_n=20):
    """
    Compute global feature importance using SHAP
    """
    if background_shap_values is None:
        return {
            "global_importance": [],
            "total_features": len(FEATURE_COLUMNS),
            "explanation": "SHAP values not available"
        }
    
    # Use mean absolute SHAP values across background dataset
    mean_abs_shap = np.mean(np.abs(background_shap_values), axis=0)
    
    # Create feature importance ranking
    importance_data = [
        {
            "feature": feature.replace("T", "U"),
            "importance": float(mean_abs_shap[i]),
            "rank": i + 1
        }
        for i, feature in enumerate(FEATURE_COLUMNS)
    ]
    
    # Sort by importance
    importance_data.sort(key=lambda x: x["importance"], reverse=True)
    
    # Add rank
    for i, item in enumerate(importance_data):
        item["rank"] = i + 1
    
    return {
        "global_importance": importance_data[:top_n],
        "total_features": len(FEATURE_COLUMNS),
        "explanation": "Features ranked by mean absolute SHAP value across all species"
    }


def shap_waterfall_data(codon, species=None):
    """
    Prepare data for waterfall plot visualization
    """
    explanation = shap_explain_codon(codon, species, detailed=False)
    
    if not explanation:
        return None
    
    # Prepare waterfall data
    features = explanation["top_features"]
    base_value = explanation["base_value"]
    
    waterfall = {
        "base_value": base_value,
        "final_prediction": explanation["prediction"],
        "steps": []
    }
    
    cumulative = base_value
    
    for feature_data in features:
        cumulative += feature_data["shap_value"]
        waterfall["steps"].append({
            "feature": feature_data["feature"],
            "shap_value": feature_data["shap_value"],
            "cumulative_value": round(cumulative, 6),
            "effect": feature_data["effect"]
        })
    
    return waterfall


def shap_force_plot_data(codon, species=None):
    """
    Prepare data for force plot visualization
    """
    explanation = shap_explain_codon(codon, species, detailed=False)
    
    if not explanation:
        return None
    
    features = explanation["top_features"]
    
    # Separate into positive and negative contributors
    positive_features = [f for f in features if f["shap_value"] > 0]
    negative_features = [f for f in features if f["shap_value"] < 0]
    
    return {
        "base_value": explanation["base_value"],
        "prediction": explanation["prediction"],
        "positive_features": positive_features,
        "negative_features": negative_features,
        "total_positive_effect": sum(f["shap_value"] for f in positive_features),
        "total_negative_effect": sum(f["shap_value"] for f in negative_features)
    }


def shap_summary_by_amino_acid(aa):
    """
    Provide SHAP summary for all codons encoding an amino acid
    """
    aa = aa.upper()
    
    if aa not in AA_TO_CODONS:
        return None
    
    codons = [c for c in AA_TO_CODONS[aa] if c in codon_df.columns]
    
    summary = {
        "amino_acid": aa,
        "codon_count": len(codons),
        "codon_explanations": []
    }
    
    for codon in codons:
        explanation = shap_explain_codon(codon, detailed=False)
        if explanation:
            summary["codon_explanations"].append({
                "codon": codon,
                "prediction": explanation["prediction"],
                "total_effect": explanation["total_effect"],
                "top_feature": explanation["top_features"][0] if explanation["top_features"] else None
            })
    
    # Rank codons by prediction value
    summary["codon_explanations"].sort(
        key=lambda x: x["prediction"], 
        reverse=True
    )
    
    return summary


# ================= ROUTES =================

@app.route("/")
def home():
    return render_template("index.html")


@app.route("/analyze", methods=["POST"])
def analyze_api():
    data = request.json

    aa = data.get("amino_acid", "")
    codon = data.get("codon", "")
    host = data.get("host_species", "")

    result, error = analyze(aa, codon)
    if error:
        return jsonify({"error": error}), 400

    response_data = {
        "codon_ranking": result["codon_ranking"],
        "selected_rank": result["selected_rank"],
        "top_species": result["top_species"],
        
        "species_specific_analysis": species_specific_analysis(aa, codon),
        "host_aware_optimization": host_aware_optimization(aa, host),
        "rare_codon_risk": rare_codon_risk(aa, host),
        "host_compatibility_score": host_compatibility_score(aa, host),

        # Enhanced SHAP explanations
        "shap_explanation": shap_explain_codon(codon, host, detailed=True) if codon else None,
        
        "model_metrics": {
            "top1_accuracy": EVAL_METRICS.get("accuracy_top1"),
            "top2_accuracy": EVAL_METRICS.get("accuracy_top2"),
            "top3_accuracy": EVAL_METRICS.get("accuracy_top3"),
            "precision": EVAL_METRICS.get("precision"),
            "recall": EVAL_METRICS.get("recall"),
            "f1_score": EVAL_METRICS.get("f1_score"),
            "loss": EVAL_METRICS.get("loss"),
            "accuracy_clean": EVAL_METRICS.get("accuracy_clean"),
            "accuracy_noisy": EVAL_METRICS.get("accuracy_noisy"),
            "accuracy_missing": EVAL_METRICS.get("accuracy_missing"),
            "accuracy_codon_only": EVAL_METRICS.get("accuracy_codon_only"),
            "accuracy_codon_bwt": EVAL_METRICS.get("accuracy_codon_bwt")
        }
    }
    
    return jsonify(response_data)


@app.route("/shap/compare", methods=["POST"])
def shap_compare_api():
    """
    Compare SHAP values for multiple codons
    """
    data = request.json
    codons = data.get("codons", [])
    species = data.get("species", None)
    
    if not codons:
        return jsonify({"error": "No codons provided"}), 400
    
    comparison = shap_compare_codons(codons, species)
    return jsonify(comparison)


@app.route("/shap/global-importance", methods=["GET"])
def shap_global_importance_api():
    """
    Get global feature importance
    """
    top_n = request.args.get("top_n", 20, type=int)
    importance = shap_global_importance(top_n)
    return jsonify(importance)


@app.route("/shap/waterfall/<codon>", methods=["GET"])
def shap_waterfall_api(codon):
    """
    Get waterfall plot data for a specific codon
    """
    species = request.args.get("species", None)
    waterfall = shap_waterfall_data(codon, species)
    
    if not waterfall:
        return jsonify({"error": "Invalid codon"}), 400
    
    return jsonify(waterfall)


@app.route("/shap/force-plot/<codon>", methods=["GET"])
def shap_force_plot_api(codon):
    """
    Get force plot data for a specific codon
    """
    species = request.args.get("species", None)
    force_data = shap_force_plot_data(codon, species)
    
    if not force_data:
        return jsonify({"error": "Invalid codon"}), 400
    
    return jsonify(force_data)


@app.route("/shap/amino-acid-summary/<aa>", methods=["GET"])
def shap_amino_acid_summary_api(aa):
    """
    Get SHAP summary for all codons of an amino acid
    """
    summary = shap_summary_by_amino_acid(aa)
    
    if not summary:
        return jsonify({"error": "Invalid amino acid"}), 400
    
    return jsonify(summary)


# ================= RUN =================

if __name__ == "__main__":
    print("\n" + "="*50)
    print("🧬 Codon Usage Tool with Explainable AI")
    print("="*50)
    print(f"✓ Server starting...")
    print(f"✓ Features: {len(FEATURE_COLUMNS)} codon features loaded")
    print(f"✓ SHAP: Ready for explanations")
    print("="*50 + "\n")
    app.run(debug=True)