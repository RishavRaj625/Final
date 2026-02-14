// ============= UTILITY FUNCTIONS =============

function showError(message) {
  const errorDiv = document.getElementById("error-message");
  errorDiv.textContent = message;
  errorDiv.style.display = "block";
  document.getElementById("loading").style.display = "none";
  document.getElementById("home-results").style.display = "none";
}

function hideError() {
  document.getElementById("error-message").style.display = "none";
}

function showLoading() {
  document.getElementById("loading").style.display = "block";
  hideError();
}

function hideLoading() {
  document.getElementById("loading").style.display = "none";
}

// ============= PAGE NAVIGATION =============

function showPage(pageName) {
  // Hide all pages
  const pages = document.querySelectorAll('.page');
  pages.forEach(page => page.classList.remove('active'));
  
  // Deactivate all nav icons
  const navIcons = document.querySelectorAll('.nav-icon');
  navIcons.forEach(icon => icon.classList.remove('active'));
  
  // Show selected page
  const selectedPage = document.getElementById(pageName);
  if (selectedPage) {
    selectedPage.classList.add('active');
  }
  
  // Activate selected nav icon
  const activeIcon = document.querySelector(`[onclick="showPage('${pageName}')"]`);
  if (activeIcon) {
    activeIcon.classList.add('active');
  }
}

// ============= IMAGE EXPANSION =============

function expandImage(img) {
  const overlay = document.getElementById("overlay");
  const imageTitle = document.getElementById("imageTitle");
  
  if (img.classList.contains("expanded")) {
    // Collapse image
    img.classList.remove("expanded");
    overlay.classList.remove("active");
    imageTitle.classList.remove("active");
    imageTitle.textContent = "";
  } else {
    // Collapse any other expanded images first
    document.querySelectorAll(".img-grid img.expanded").forEach(i => {
      i.classList.remove("expanded");
    });
    
    // Expand this image
    img.classList.add("expanded");
    overlay.classList.add("active");
    
    // Set image title
    const altText = img.getAttribute("alt") || "Image";
    imageTitle.textContent = altText;
    imageTitle.classList.add("active");
  }
}

// Close expanded image when clicking overlay
document.addEventListener('DOMContentLoaded', function() {
  const overlay = document.getElementById("overlay");
  if (overlay) {
    overlay.addEventListener('click', function() {
      document.querySelectorAll(".img-grid img.expanded").forEach(img => {
        img.classList.remove("expanded");
      });
      overlay.classList.remove("active");
      const imageTitle = document.getElementById("imageTitle");
      if (imageTitle) {
        imageTitle.classList.remove("active");
        imageTitle.textContent = "";
      }
    });
  }
});

// ============= MAIN ANALYSIS FUNCTION =============

async function analyze() {
  // Get input values
  const aa = document.getElementById("aa").value.trim().toUpperCase();
  const codon = document.getElementById("codon").value.trim().toUpperCase();
  const host = document.getElementById("host").value.trim();
  
  // Validate input
  if (!aa) {
    showError("Please enter an amino acid (single letter code)");
    return;
  }
  
  // Show loading
  showLoading();
  hideError();
  
  try {
    // Make API call
    const response = await fetch("/analyze", {
      method: "POST",
      headers: {
        "Content-Type": "application/json"
      },
      body: JSON.stringify({
        amino_acid: aa,
        codon: codon,
        host_species: host
      })
    });
    
    // Check if response is ok
    if (!response.ok) {
      const errorData = await response.json();
      throw new Error(errorData.error || "Server error occurred");
    }
    
    const data = await response.json();
    
    // Hide loading and show results
    hideLoading();
    displayResults(data);
    
    // Update metrics page if model_metrics exists
    if (data.model_metrics) {
      updateMetricsPage(data.model_metrics);
    }
    
  } catch (error) {
    console.error("Analysis error:", error);
    showError(error.message || "Failed to analyze. Please check your inputs and try again.");
  }
}

// ============= DISPLAY RESULTS =============

function displayResults(data) {
  // Show results container
  document.getElementById("home-results").style.display = "block";
  
  // 1. Codon Ranking Table
  displayCodonRanking(data.codon_ranking, data.selected_rank);
  
  // 2. Top Species
  displayTopSpecies(data.top_species);
  
  // 3. Species Preferences
  if (data.species_specific_analysis) {
    displaySpeciesPreferences(data.species_specific_analysis);
  }
  
  // 4. Host Optimization
  if (data.host_aware_optimization) {
    displayHostOptimization(data.host_aware_optimization);
  }
  
  // 5. Codon Bias
  if (data.codon_bias_score) {
    displayCodonBias(data.codon_bias_score);
  }
  
  // 6. Kingdom Comparison
  if (data.cross_kingdom_comparison) {
    displayKingdomComparison(data.cross_kingdom_comparison);
  }
  
  // 7. SHAP Explanation
  if (data.shap_explanation) {
    displayShapExplanation(data.shap_explanation);
  }
}

// ============= INDIVIDUAL DISPLAY FUNCTIONS =============

function displayCodonRanking(ranking, selectedRank) {
  const tbody = document.querySelector("#rankingTable tbody");
  tbody.innerHTML = "";
  
  if (!ranking || ranking.length === 0) {
    tbody.innerHTML = '<tr><td colspan="4">No data available</td></tr>';
    return;
  }
  
  ranking.forEach(row => {
    const tr = document.createElement("tr");
    
    // Highlight if: (1) this is rank 1, OR (2) this is the selected codon
    const shouldHighlight = row.rank === 1 || (selectedRank && row.rank === selectedRank);
    if (shouldHighlight) {
      tr.classList.add("highlight");
    }
    
    tr.innerHTML = `
      <td>${row.rank}</td>
      <td><strong>${row.codon}</strong></td>
      <td>${row.frequency.toFixed(4)}</td>
      <td>${row.ml_weight.toFixed(4)}</td>
    `;
    tbody.appendChild(tr);
  });
}

function displayTopSpecies(species) {
  const ul = document.getElementById("species");
  ul.innerHTML = "";
  
  if (!species || species.length === 0) {
    ul.innerHTML = '<li>No species data available</li>';
    return;
  }
  
  species.forEach((s, index) => {
    const li = document.createElement("li");
    
    // Add special class for top species
    if (index === 0) {
      li.classList.add("top-species");
    }
    
    li.textContent = `${s.SPECIESNAME || s.species} (Score: ${(s.SCORE || s.usage).toFixed(4)})`;
    ul.appendChild(li);
  });
}

function displaySpeciesPreferences(analysis) {
  if (!analysis) return;
  
  // High preference species
  const highList = document.getElementById("speciesHigh");
  highList.innerHTML = "";
  
  if (analysis.top_species && analysis.top_species.length > 0) {
    analysis.top_species.forEach(s => {
      const li = document.createElement("li");
      li.innerHTML = `
        <span>${s.SPECIESNAME}</span>
        <span class="badge">${s.PREFERENCE_SCORE.toFixed(3)}</span>
      `;
      highList.appendChild(li);
    });
  } else {
    highList.innerHTML = '<li>No data available</li>';
  }
  
  // Low preference species
  const lowList = document.getElementById("speciesLow");
  lowList.innerHTML = "";
  
  if (analysis.bottom_species && analysis.bottom_species.length > 0) {
    analysis.bottom_species.forEach(s => {
      const li = document.createElement("li");
      li.innerHTML = `
        <span>${s.SPECIESNAME}</span>
        <span class="badge">${s.PREFERENCE_SCORE.toFixed(3)}</span>
      `;
      lowList.appendChild(li);
    });
  } else {
    lowList.innerHTML = '<li>No data available</li>';
  }
  
  // Explanation
  const explain = document.getElementById("speciesExplain");
  explain.textContent = analysis.explanation || "Species with high preference frequently use this codon, while low preference species rarely use it.";
}

function displayHostOptimization(optimization) {
  const hostResult = document.getElementById("hostResult");
  const container = document.getElementById("hostRankingContainer");
  const tbody = document.querySelector("#hostRankingTable tbody");
  
  if (!optimization) {
    hostResult.innerHTML = `
      <p style="color:#888;">💡 <strong>Enter a host species name</strong> above to see optimized codon usage.</p>
    `;
    container.style.display = "none";
    return;
  }
  
  // Check if data was found
  if (!optimization.found) {
    hostResult.innerHTML = `
      <div style="padding:15px;background:rgba(255,165,0,0.1);border-left:4px solid #FFA500;border-radius:5px;">
        <strong style="color:#FFA500;">ℹ️ ${optimization.message || 'No host specified'}</strong><br>
        <em style="color:#ccc;font-size:13px;">Examples: "Escherichia coli", "E. coli", "coli", "Human", "Yeast", "Bacillus subtilis"</em>
      </div>
    `;
    container.style.display = "none";
    return;
  }
  
  // Data found - show results
  hostResult.innerHTML = `
    <div style="padding:15px;background:rgba(76,175,80,0.1);border-left:4px solid #4CAF50;border-radius:5px;">
      <strong style="color:#4CAF50;">✓ Host Species Found:</strong> ${optimization.host_species}<br>
      <strong style="color:#d4af37;">🎯 Optimal Codon:</strong> <span style="color:#4CAF50;font-size:20px;font-weight:bold;">${optimization.optimal_codon}</span><br>
      <em style="color:#ccc;font-size:13px;">This is the most frequently used codon for this amino acid in ${optimization.host_species}</em>
    </div>
  `;
  
  if (optimization.codon_ranking && optimization.codon_ranking.length > 0) {
    container.style.display = "block";
    tbody.innerHTML = "";
    
    optimization.codon_ranking.forEach((item, index) => {
      const tr = document.createElement("tr");
      
      // Highlight rank 1
      if (index === 0) {
        tr.classList.add("highlight");
      }
      
      tr.innerHTML = `
        <td>${index + 1}</td>
        <td><strong>${item[0]}</strong></td>
        <td>${item[1].toFixed(4)}</td>
      `;
      tbody.appendChild(tr);
    });
  } else {
    container.style.display = "none";
  }
}

function displayCodonBias(biasData) {
  const biasResult = document.getElementById("biasResult");
  const container = document.getElementById("biasTableContainer");
  const tbody = document.querySelector("#biasTable tbody");
  
  if (!biasData) {
    biasResult.innerHTML = `
      <p style="color:#888;">💡 <strong>Enter a specific codon</strong> above to see bias analysis.</p>
    `;
    container.style.display = "none";
    return;
  }
  
  // Check if data was found
  if (!biasData.found) {
    biasResult.innerHTML = `
      <div style="padding:15px;background:rgba(255,165,0,0.1);border-left:4px solid #FFA500;border-radius:5px;">
        <strong style="color:#FFA500;">ℹ️ ${biasData.message || 'No codon specified'}</strong><br>
        <em style="color:#ccc;font-size:13px;">Enter a codon in RNA format (e.g., UUA, GCC, UAA) in the "Codon" field above</em>
      </div>
    `;
    container.style.display = "none";
    return;
  }
  
  // Data found - show results
  biasResult.innerHTML = `
    <div style="padding:15px;background:rgba(76,175,80,0.1);border-left:4px solid #4CAF50;border-radius:5px;">
      <strong style="color:#4CAF50;">✓ Analyzing Codon:</strong> <span style="font-size:18px;font-weight:bold;color:#d4af37;">${biasData.codon}</span><br>
      <strong>Global Average Usage:</strong> ${biasData.global_average.toFixed(4)}<br>
      <em style="color:#ccc;font-size:13px;">Bias score = (Species usage) / (Global average). Higher values indicate stronger preference.</em>
    </div>
  `;
  
  if (biasData.top_bias_species && biasData.top_bias_species.length > 0) {
    container.style.display = "block";
    tbody.innerHTML = "";
    
    biasData.top_bias_species.forEach((s, index) => {
      const tr = document.createElement("tr");
      
      // Highlight top species
      if (index === 0) {
        tr.classList.add("highlight");
      }
      
      tr.innerHTML = `
        <td>${s.SPECIESNAME}</td>
        <td><strong>${s.bias.toFixed(2)}x</strong></td>
      `;
      tbody.appendChild(tr);
    });
  } else {
    container.style.display = "none";
  }
}

function displayKingdomComparison(kingdoms) {
  const tbody = document.querySelector("#kingdomTable tbody");
  tbody.innerHTML = "";
  
  if (!kingdoms || kingdoms.length === 0) {
    tbody.innerHTML = '<tr><td colspan="2">No kingdom data available</td></tr>';
    return;
  }
  
  // Sort by usage descending
  kingdoms.sort((a, b) => {
    const aVal = a[Object.keys(a)[1]] || 0;
    const bVal = b[Object.keys(b)[1]] || 0;
    return bVal - aVal;
  });
  
  kingdoms.forEach((k, index) => {
    const tr = document.createElement("tr");
    
    // Highlight highest usage kingdom
    if (index === 0) {
      tr.classList.add("highlight");
    }
    
    const kingdomName = k.KINGDOM || k.kingdom;
    const usage = k[Object.keys(k)[1]];
    
    tr.innerHTML = `
      <td>${kingdomName}</td>
      <td>${usage.toFixed(4)}</td>
    `;
    tbody.appendChild(tr);
  });
}

function displayShapExplanation(shap) {
  const chartDiv = document.getElementById("shapChart");
  const tbody = document.querySelector("#shapTable tbody");
  
  if (!chartDiv || !tbody) return;
  
  chartDiv.innerHTML = "";
  tbody.innerHTML = "";
  
  if (!shap || shap.length === 0) {
    chartDiv.innerHTML = '<p style="text-align:center;color:#888;">No SHAP data available</p>';
    tbody.innerHTML = '<tr><td colspan="3">No data available</td></tr>';
    return;
  }
  
  // Sort by absolute impact (highest influence first)
  const sortedShap = [...shap].sort((a, b) => Math.abs(b.impact) - Math.abs(a.impact));
  
  // Create bar chart
  sortedShap.forEach((item, index) => {
    const barContainer = document.createElement("div");
    barContainer.style.margin = "12px 0";
    barContainer.style.display = "flex";
    barContainer.style.alignItems = "center";
    barContainer.style.gap = "12px";
    
    const label = document.createElement("strong");
    label.textContent = item.feature || item.codon || item.name || "Unknown";
    label.style.minWidth = "60px";
    label.style.color = "#d4af37";
    label.style.fontSize = "16px";
    
    const barWrapper = document.createElement("div");
    barWrapper.style.flex = "1";
    barWrapper.style.height = "30px";
    barWrapper.style.background = "#1a1a1a";
    barWrapper.style.position = "relative";
    barWrapper.style.borderRadius = "5px";
    barWrapper.style.overflow = "hidden";
    barWrapper.style.border = "1px solid #333";
    
    const impact = item.impact || item.value || 0;
    const absImpact = Math.abs(impact);
    const maxImpact = Math.max(...sortedShap.map(s => Math.abs(s.impact || s.value || 0)));
    const barWidth = maxImpact > 0 ? (absImpact / maxImpact * 100) : 0;
    
    const bar = document.createElement("div");
    bar.style.width = barWidth + "%";
    bar.style.height = "100%";
    bar.style.background = impact > 0 ? 
      "linear-gradient(90deg, #4CAF50, #66BB6A)" : 
      "linear-gradient(90deg, #FF5252, #FF7961)";
    bar.style.transition = "width 0.5s ease";
    bar.style.boxShadow = impact > 0 ? 
      "0 0 10px rgba(76, 175, 80, 0.5)" : 
      "0 0 10px rgba(255, 82, 82, 0.5)";
    
    const value = document.createElement("span");
    value.textContent = `${impact > 0 ? '+' : ''}${impact.toFixed(4)}`;
    value.style.marginLeft = "10px";
    value.style.color = impact > 0 ? "#4CAF50" : "#FF5252";
    value.style.fontWeight = "bold";
    value.style.minWidth = "90px";
    value.style.fontSize = "14px";
    
    barWrapper.appendChild(bar);
    barContainer.appendChild(label);
    barContainer.appendChild(barWrapper);
    barContainer.appendChild(value);
    chartDiv.appendChild(barContainer);
    
    // Add to table with interpretation
    const tr = document.createElement("tr");
    
    // Highlight top influence
    if (index === 0) {
      tr.classList.add("highlight");
    }
    
    let interpretation = "";
    if (impact > 0) {
      if (absImpact > 0.02) {
        interpretation = "Strongly recommended";
      } else if (absImpact > 0.01) {
        interpretation = "Moderately recommended";
      } else {
        interpretation = "Slightly favored";
      }
    } else {
      if (absImpact > 0.02) {
        interpretation = "Strongly discouraged";
      } else if (absImpact > 0.01) {
        interpretation = "Moderately discouraged";
      } else {
        interpretation = "Slightly disfavored";
      }
    }
    
    tr.innerHTML = `
      <td><strong>${item.feature || item.codon || "Unknown"}</strong></td>
      <td style="color: ${impact > 0 ? '#4CAF50' : '#FF5252'}; font-weight: bold;">
        ${impact > 0 ? '+' : ''}${impact.toFixed(4)}
      </td>
      <td style="color: #ccc;">${interpretation}</td>
    `;
    tbody.appendChild(tr);
  });
}

// ============= UPDATE METRICS PAGE =============

function updateMetricsPage(metrics) {
  if (!metrics) return;
  
  console.log("Updating metrics with:", metrics);
  
  // Accuracy Comparison
  const accCodon = document.getElementById("accCodon");
  const accBWT = document.getElementById("accBWT");
  
  if (accCodon && metrics.accuracy_codon_only !== undefined) {
    accCodon.textContent = (metrics.accuracy_codon_only * 100).toFixed(2) + "%";
  }
  if (accBWT && (metrics.accuracy_codon_bwt !== undefined || metrics.accuracy_codon_BWT !== undefined)) {
    const val = metrics.accuracy_codon_bwt || metrics.accuracy_codon_BWT;
    accBWT.textContent = (val * 100).toFixed(2) + "%";
  }
  
  // Model Evaluation Metrics - handle both naming conventions
  const metricsTableBody = document.querySelector("#metricsTable tbody");
  if (metricsTableBody) {
    const rows = metricsTableBody.querySelectorAll("tr");
    
    // Support both top1_accuracy and accuracy_top1 naming
    const top1 = metrics.top1_accuracy || metrics.accuracy_top1;
    const top2 = metrics.top2_accuracy || metrics.accuracy_top2;
    const top3 = metrics.top3_accuracy || metrics.accuracy_top3;
    
    const metricsMap = {
      0: top1 !== undefined ? (top1 * 100).toFixed(2) + "%" : "–",
      1: top2 !== undefined ? (top2 * 100).toFixed(2) + "%" : "–",
      2: top3 !== undefined ? (top3 * 100).toFixed(2) + "%" : "–",
      3: metrics.precision !== undefined ? (metrics.precision * 100).toFixed(2) + "%" : "–",
      4: metrics.recall !== undefined ? (metrics.recall * 100).toFixed(2) + "%" : "–",
      5: metrics.f1_score !== undefined ? metrics.f1_score.toFixed(4) : "–",
      6: metrics.loss !== undefined ? metrics.loss.toFixed(4) : "–",
      7: top1 !== undefined ? ((1 - top1) * 100).toFixed(2) + "%" : "–"
    };
    
    rows.forEach((row, index) => {
      const cells = row.querySelectorAll("td");
      if (cells.length > 1 && metricsMap[index] !== undefined) {
        cells[1].textContent = metricsMap[index];
      }
    });
  }
  
  // Robustness Evaluation
  const accClean = document.getElementById("accClean");
  const accNoisy = document.getElementById("accNoisy");
  const accMissing = document.getElementById("accMissing");
  
  if (accClean && metrics.accuracy_clean !== undefined) {
    accClean.textContent = (metrics.accuracy_clean * 100).toFixed(2) + "%";
  }
  if (accNoisy && metrics.accuracy_noisy !== undefined) {
    accNoisy.textContent = (metrics.accuracy_noisy * 100).toFixed(2) + "%";
  }
  if (accMissing && metrics.accuracy_missing !== undefined) {
    accMissing.textContent = (metrics.accuracy_missing * 100).toFixed(2) + "%";
  }
}

// ============= KEYBOARD SUPPORT =============

document.addEventListener('DOMContentLoaded', function() {
  // Allow Enter key to trigger analysis
  const inputs = ['aa', 'codon', 'host'];
  inputs.forEach(id => {
    const input = document.getElementById(id);
    if (input) {
      input.addEventListener('keypress', function(e) {
        if (e.key === 'Enter') {
          analyze();
        }
      });
    }
  });
});

// ============= NETWORK BACKGROUND =============

const canvas = document.getElementById('network-bg');
if (canvas) {
  const ctx = canvas.getContext('2d');

  canvas.width = window.innerWidth;
  canvas.height = window.innerHeight;

  window.addEventListener('resize', function() {
    canvas.width = window.innerWidth;
    canvas.height = window.innerHeight;
  });

  // Simple particle animation
  const particles = [];
  const particleCount = 50;

  for (let i = 0; i < particleCount; i++) {
    particles.push({
      x: Math.random() * canvas.width,
      y: Math.random() * canvas.height,
      vx: (Math.random() - 0.5) * 0.5,
      vy: (Math.random() - 0.5) * 0.5,
      radius: Math.random() * 2 + 1
    });
  }

  function animate() {
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    
    particles.forEach(p => {
      p.x += p.vx;
      p.y += p.vy;
      
      if (p.x < 0 || p.x > canvas.width) p.vx *= -1;
      if (p.y < 0 || p.y > canvas.height) p.vy *= -1;
      
      ctx.beginPath();
      ctx.arc(p.x, p.y, p.radius, 0, Math.PI * 2);
      ctx.fillStyle = 'rgba(212, 175, 55, 0.5)';
      ctx.fill();
    });
    
    requestAnimationFrame(animate);
  }

  animate();
}