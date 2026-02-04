// Network Background Animation
const canvas = document.getElementById('network-bg');
const ctx = canvas.getContext('2d');
canvas.width = window.innerWidth;
canvas.height = window.innerHeight;

const particles = [];
const particleCount = window.innerWidth < 768 ? 40 : 80;

class Particle {
  constructor() {
    this.x = Math.random() * canvas.width;
    this.y = Math.random() * canvas.height;
    this.vx = (Math.random() - 0.5) * 0.5;
    this.vy = (Math.random() - 0.5) * 0.5;
    this.radius = Math.random() * 2 + 1;
  }

  update() {
    this.x += this.vx;
    this.y += this.vy;

    if (this.x < 0 || this.x > canvas.width) this.vx *= -1;
    if (this.y < 0 || this.y > canvas.height) this.vy *= -1;
  }

  draw() {
    ctx.beginPath();
    ctx.arc(this.x, this.y, this.radius, 0, Math.PI * 2);
    ctx.fillStyle = 'rgba(212, 175, 55, 0.5)';
    ctx.fill();
  }
}

for (let i = 0; i < particleCount; i++) {
  particles.push(new Particle());
}

function animate() {
  ctx.clearRect(0, 0, canvas.width, canvas.height);

  particles.forEach(p => {
    p.update();
    p.draw();
  });

  // Draw connections
  particles.forEach((p1, i) => {
    particles.slice(i + 1).forEach(p2 => {
      const dx = p1.x - p2.x;
      const dy = p1.y - p2.y;
      const dist = Math.sqrt(dx * dx + dy * dy);

      if (dist < 150) {
        ctx.beginPath();
        ctx.moveTo(p1.x, p1.y);
        ctx.lineTo(p2.x, p2.y);
        ctx.strokeStyle = `rgba(212, 175, 55, ${0.2 * (1 - dist / 150)})`;
        ctx.lineWidth = 1;
        ctx.stroke();
      }
    });
  });

  requestAnimationFrame(animate);
}

animate();

window.addEventListener('resize', () => {
  canvas.width = window.innerWidth;
  canvas.height = window.innerHeight;
});

// Global variable to store current SHAP data
let currentShapData = null;
let currentCodon = null;
let currentHost = null;

// Chart instances (for cleanup)
let shapBarChartInstance = null;
let waterfallChartInstance = null;
let globalImportanceChartInstance = null;

// Navigation
function showPage(pageId) {
  document.querySelectorAll('.page').forEach(p => p.classList.remove('active'));
  document.querySelectorAll('.nav-icon').forEach(n => n.classList.remove('active'));
  
  document.getElementById(pageId).classList.add('active');
  event.target.closest('.nav-icon').classList.add('active');
  
  // If switching to SHAP page, display current SHAP data if available
  if (pageId === 'shap' && currentShapData) {
    displayShapExplanation(currentShapData);
  }
}

// Image Expansion
let expandedImage = null;
const overlay = document.getElementById('overlay');
const imageTitle = document.getElementById('imageTitle');

function expandImage(img) {
  if (expandedImage && expandedImage !== img) {
    expandedImage.classList.remove('expanded');
  }

  if (img.classList.contains('expanded')) {
    img.classList.remove('expanded');
    overlay.classList.remove('active');
    imageTitle.classList.remove('active');
    expandedImage = null;
  } else {
    img.classList.add('expanded');
    overlay.classList.add('active');
    imageTitle.textContent = img.alt;
    imageTitle.classList.add('active');
    expandedImage = img;
  }
}

overlay.addEventListener('click', () => {
  if (expandedImage) {
    expandedImage.classList.remove('expanded');
    overlay.classList.remove('active');
    imageTitle.classList.remove('active');
    expandedImage = null;
  }
});

// Analyze Function
async function analyze() {
  const aa = document.getElementById("aa").value.trim().toUpperCase();
  const codon = document.getElementById("codon").value.trim().toUpperCase();
  const host = document.getElementById("host").value.trim();

  // Store current inputs
  currentCodon = codon;
  currentHost = host;

  // Show loading
  document.getElementById('loading').style.display = 'block';
  document.getElementById('error-message').style.display = 'none';
  document.getElementById('home-results').style.display = 'none';

  try {
    const res = await fetch("/analyze", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ amino_acid: aa, codon: codon, host_species: host })
    });

    if (!res.ok) {
      throw new Error(`Server error: ${res.status} ${res.statusText}`);
    }

    const data = await res.json();
    
    // Hide loading
    document.getElementById('loading').style.display = 'none';
    
    if (data.error) {
      document.getElementById('error-message').textContent = data.error;
      document.getElementById('error-message').style.display = 'block';
      return;
    }

    // Store SHAP data globally
    currentShapData = data.shap_explanation;

    // Show results
    document.getElementById('home-results').style.display = 'block';

    // Codon ranking
    let rHTML = '';
    data.codon_ranking.forEach(r => {
      rHTML += `<tr class="${r.codon === codon ? 'highlight' : ''}">
        <td>${r.rank}</td><td>${r.codon}</td>
        <td>${r.frequency.toFixed(4)}</td>
        <td>${r.ml_weight.toFixed(2)}</td></tr>`;
    });
    document.querySelector("#rankingTable tbody").innerHTML = rHTML;

    // Top species
    document.getElementById("species").innerHTML =
      data.top_species.map((s, index) =>
        `<li class="${index === 0 ? 'top-species' : ''}">${s.SPECIESNAME} — ${s.SCORE.toFixed(4)}</li>`
      ).join("");

    // Species-specific analysis
    if (data.species_specific_analysis) {
      document.getElementById("speciesHigh").innerHTML =
        data.species_specific_analysis.top_species.map(s =>
          `<li>${s.SPECIESNAME} <span class="badge">${(s.PREFERENCE_SCORE * 100).toFixed(1)}%</span></li>`
        ).join("");

      document.getElementById("speciesLow").innerHTML =
        data.species_specific_analysis.bottom_species.map(s =>
          `<li>${s.SPECIESNAME} <span class="badge">${(s.PREFERENCE_SCORE * 100).toFixed(1)}%</span></li>`
        ).join("");

      document.getElementById("speciesExplain").innerHTML =
        `<strong>Analysis:</strong> ${data.species_specific_analysis.explanation}`;
    }

    // Host-aware optimization
    const hostResultDiv = document.getElementById("hostResult");
    const hostRankingContainer = document.getElementById('hostRankingContainer');
    
    if (data.host_aware_optimization) {
      const hostOpt = data.host_aware_optimization;
      hostResultDiv.innerHTML =
        `<strong>Optimal codon for ${hostOpt.host_species}:</strong> 
         <span style="color: #d4af37; font-size: 16px; font-weight: bold;">${hostOpt.optimal_codon}</span>
         <br><br>This codon shows the highest usage frequency in the selected host organism.`;

      if (hostOpt.codon_ranking && hostOpt.codon_ranking.length > 0) {
        hostRankingContainer.style.display = 'block';
        let hostRankHTML = '';
        hostOpt.codon_ranking.forEach(([codon, usage], index) => {
          hostRankHTML += `<tr class="${index === 0 ? 'highlight' : ''}">
            <td>${index + 1}</td>
            <td>${codon}</td>
            <td>${usage.toFixed(4)}</td>
          </tr>`;
        });
        document.querySelector("#hostRankingTable tbody").innerHTML = hostRankHTML;
      } else {
        hostRankingContainer.style.display = 'none';
      }
    } else {
      hostResultDiv.innerHTML =
        `<strong>No host species selected.</strong><br>
         Enter a host species name (e.g., "Escherichia coli") in the input field above to see optimized codon recommendations for that organism.`;
      hostRankingContainer.style.display = 'none';
    }

    // Metrics
    document.getElementById("accCodon").innerText =
      (data.model_metrics.accuracy_codon_only * 100).toFixed(2) + "%";
    document.getElementById("accBWT").innerText =
      (data.model_metrics.accuracy_codon_bwt * 100).toFixed(2) + "%";

    const m = data.model_metrics;
    const rows = document.querySelectorAll("#metricsTable tbody tr td:nth-child(2)");
    rows[0].innerText = m.top1_accuracy.toFixed(4);
    rows[1].innerText = (m.top2_accuracy ?? 0).toFixed(4);
    rows[2].innerText = m.top3_accuracy.toFixed(4);
    rows[3].innerText = m.precision.toFixed(4);
    rows[4].innerText = m.recall.toFixed(4);
    rows[5].innerText = m.f1_score.toFixed(4);
    rows[6].innerText = m.loss.toFixed(4);
    rows[7].innerText = (1 - m.top1_accuracy).toFixed(4);

    document.getElementById("accClean").innerText =
      (m.accuracy_clean * 100).toFixed(2) + "%";
    document.getElementById("accNoisy").innerText =
      (m.accuracy_noisy * 100).toFixed(2) + "%";
    document.getElementById("accMissing").innerText =
      (m.accuracy_missing * 100).toFixed(2) + "%";
    
  } catch (error) {
    document.getElementById('loading').style.display = 'none';
    document.getElementById('error-message').textContent = 
      `Error: ${error.message}. Make sure Flask server is running on the correct port.`;
    document.getElementById('error-message').style.display = 'block';
    console.error('Fetch error:', error);
  }
}

// Display SHAP Explanation
function displayShapExplanation(shapData) {
  if (!shapData) {
    document.getElementById('shap-results').style.display = 'none';
    document.getElementById('no-shap-data').style.display = 'block';
    return;
  }

  document.getElementById('shap-results').style.display = 'block';
  document.getElementById('no-shap-data').style.display = 'none';

  // Summary
  const summary = `
    <strong>Codon:</strong> ${shapData.codon}<br>
    <strong>Base Value (Average):</strong> ${shapData.base_value.toFixed(4)}<br>
    <strong>Prediction:</strong> ${shapData.prediction.toFixed(4)}<br>
    <strong>Total Effect:</strong> ${shapData.total_effect > 0 ? '+' : ''}${shapData.total_effect.toFixed(4)}
    ${currentHost ? `<br><strong>Context:</strong> ${currentHost}` : '<br><strong>Context:</strong> Global'}
  `;
  document.getElementById('shap-summary').innerHTML = summary;

  // Top Features Table
  let featuresHTML = '';
  shapData.top_features.forEach((f, index) => {
    const effectClass = f.effect === 'positive' ? 'positive-effect' : 'negative-effect';
    featuresHTML += `<tr>
      <td>${index + 1}</td>
      <td><strong>${f.feature}</strong></td>
      <td class="${effectClass}">${f.shap_value.toFixed(6)}</td>
      <td><span class="badge ${effectClass}">${f.effect}</span></td>
      <td>${f.percentage}%</td>
    </tr>`;
  });
  document.querySelector("#shapFeaturesTable tbody").innerHTML = featuresHTML;

  // Create SHAP Bar Chart
  createShapBarChart(shapData.top_features);

  // Create Waterfall Chart
  createWaterfallChart(shapData);

  // Display Interactions if available
  if (shapData.interactions && shapData.interactions.length > 0) {
    document.getElementById('interactionsSection').style.display = 'block';
    let interactionsHTML = '';
    shapData.interactions.forEach(interaction => {
      interactionsHTML += `<tr>
        <td>${interaction.feature}</td>
        <td>${interaction.interaction_value.toFixed(6)}</td>
      </tr>`;
    });
    document.querySelector("#interactionsTable tbody").innerHTML = interactionsHTML;
  } else {
    document.getElementById('interactionsSection').style.display = 'none';
  }

  // Display Interpretation
  if (shapData.interpretation) {
    document.getElementById('shap-interpretation').innerHTML = shapData.interpretation;
  }
}

// Create SHAP Bar Chart
function createShapBarChart(features) {
  const ctx = document.getElementById('shapBarChart');
  
  // Destroy previous chart if exists
  if (shapBarChartInstance) {
    shapBarChartInstance.destroy();
  }

  const labels = features.map(f => f.feature);
  const values = features.map(f => f.shap_value);
  const colors = values.map(v => v > 0 ? 'rgba(75, 192, 75, 0.7)' : 'rgba(255, 99, 99, 0.7)');

  shapBarChartInstance = new Chart(ctx, {
    type: 'bar',
    data: {
      labels: labels,
      datasets: [{
        label: 'SHAP Value',
        data: values,
        backgroundColor: colors,
        borderColor: colors.map(c => c.replace('0.7', '1')),
        borderWidth: 2
      }]
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      indexAxis: 'y',
      plugins: {
        legend: {
          display: false
        },
        title: {
          display: true,
          text: 'Feature Contribution to Prediction',
          color: '#d4af37',
          font: {
            size: 16
          }
        },
        tooltip: {
          callbacks: {
            label: function(context) {
              return `SHAP Value: ${context.parsed.x.toFixed(6)}`;
            }
          }
        }
      },
      scales: {
        x: {
          title: {
            display: true,
            text: 'SHAP Value',
            color: '#d4af37'
          },
          grid: {
            color: 'rgba(212, 175, 55, 0.1)'
          },
          ticks: {
            color: '#ccc'
          }
        },
        y: {
          grid: {
            color: 'rgba(212, 175, 55, 0.1)'
          },
          ticks: {
            color: '#ccc',
            font: {
              size: 11
            }
          }
        }
      }
    }
  });
}

// Create Waterfall Chart
function createWaterfallChart(shapData) {
  const ctx = document.getElementById('waterfallChart');
  
  // Destroy previous chart if exists
  if (waterfallChartInstance) {
    waterfallChartInstance.destroy();
  }

  // Prepare waterfall data
  const baseValue = shapData.base_value;
  const features = shapData.top_features.slice(0, 8); // Top 8 features
  
  let cumulative = baseValue;
  const labels = ['Base Value'];
  const data = [baseValue];
  const colors = ['rgba(100, 100, 100, 0.7)'];

  features.forEach(f => {
    labels.push(f.feature);
    cumulative += f.shap_value;
    data.push(cumulative);
    colors.push(f.shap_value > 0 ? 'rgba(75, 192, 75, 0.7)' : 'rgba(255, 99, 99, 0.7)');
  });

  labels.push('Final');
  data.push(shapData.prediction);
  colors.push('rgba(212, 175, 55, 0.7)');

  waterfallChartInstance = new Chart(ctx, {
    type: 'line',
    data: {
      labels: labels,
      datasets: [{
        label: 'Cumulative Value',
        data: data,
        borderColor: '#d4af37',
        backgroundColor: colors,
        fill: false,
        tension: 0.1,
        pointRadius: 6,
        pointBackgroundColor: colors,
        pointBorderColor: '#fff',
        pointBorderWidth: 2
      }]
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      plugins: {
        legend: {
          display: false
        },
        title: {
          display: true,
          text: 'How Features Build to Final Prediction',
          color: '#d4af37',
          font: {
            size: 16
          }
        },
        tooltip: {
          callbacks: {
            label: function(context) {
              return `Value: ${context.parsed.y.toFixed(6)}`;
            }
          }
        }
      },
      scales: {
        x: {
          grid: {
            color: 'rgba(212, 175, 55, 0.1)'
          },
          ticks: {
            color: '#ccc',
            maxRotation: 45,
            minRotation: 45,
            font: {
              size: 10
            }
          }
        },
        y: {
          title: {
            display: true,
            text: 'Prediction Value',
            color: '#d4af37'
          },
          grid: {
            color: 'rgba(212, 175, 55, 0.1)'
          },
          ticks: {
            color: '#ccc'
          }
        }
      }
    }
  });
}

// Load Global Importance
async function loadGlobalImportance() {
  try {
    const res = await fetch("/shap/global-importance?top_n=20");
    const data = await res.json();

    if (data.global_importance) {
      document.getElementById('global-importance-container').style.display = 'block';
      createGlobalImportanceChart(data.global_importance);
    }
  } catch (error) {
    console.error('Error loading global importance:', error);
    alert('Failed to load global importance data.');
  }
}

// Create Global Importance Chart
function createGlobalImportanceChart(importanceData) {
  const ctx = document.getElementById('globalImportanceChart');
  
  // Destroy previous chart if exists
  if (globalImportanceChartInstance) {
    globalImportanceChartInstance.destroy();
  }

  const labels = importanceData.map(d => d.feature);
  const values = importanceData.map(d => d.importance);

  globalImportanceChartInstance = new Chart(ctx, {
    type: 'bar',
    data: {
      labels: labels,
      datasets: [{
        label: 'Mean Absolute SHAP Value',
        data: values,
        backgroundColor: 'rgba(212, 175, 55, 0.7)',
        borderColor: 'rgba(212, 175, 55, 1)',
        borderWidth: 2
      }]
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      indexAxis: 'y',
      plugins: {
        legend: {
          display: false
        },
        title: {
          display: true,
          text: 'Top 20 Most Important Features (Global)',
          color: '#d4af37',
          font: {
            size: 16
          }
        },
        tooltip: {
          callbacks: {
            label: function(context) {
              return `Importance: ${context.parsed.x.toFixed(6)}`;
            }
          }
        }
      },
      scales: {
        x: {
          title: {
            display: true,
            text: 'Mean Absolute SHAP Value',
            color: '#d4af37'
          },
          grid: {
            color: 'rgba(212, 175, 55, 0.1)'
          },
          ticks: {
            color: '#ccc'
          }
        },
        y: {
          grid: {
            color: 'rgba(212, 175, 55, 0.1)'
          },
          ticks: {
            color: '#ccc',
            font: {
              size: 10
            }
          }
        }
      }
    }
  });
}

// Test server connection on page load
window.addEventListener('DOMContentLoaded', async () => {
  console.log('Page loaded. Testing server connection...');
  
  try {
    const testRes = await fetch('/analyze', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ amino_acid: 'L', codon: 'UUA' })
    });
    
    if (testRes.ok) {
      console.log('✓ Server connection successful');
    } else {
      console.warn('⚠ Server returned:', testRes.status, testRes.statusText);
    }
  } catch (err) {
    console.error('✗ Cannot connect to server:', err.message);
    console.error('Make sure Flask is running on http://localhost:5000 or update the fetch URL');
  }
});