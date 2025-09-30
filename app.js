// Global variables
let Module = null;
let currentData = null;
let analysisResults = null;
let componentChart = null;
let fitChart = null;
let rawDataChart = null;
let componentTraceChart = null;
let logPosteriorChart = null;
let isRunning = false;
let currentIteration = 0;
let totalIterations = 0;
let mcmcState = null;

// Initialize the application
document.addEventListener('DOMContentLoaded', function() {
    initializeFileUpload();
    initializeWebAssembly();
});

// Initialize file upload functionality
function initializeFileUpload() {
    const fileUpload = document.getElementById('fileUpload');
    const fileInput = document.getElementById('fileInput');
    const fileInfo = document.getElementById('fileInfo');

    fileUpload.addEventListener('click', () => fileInput.click());
    fileUpload.addEventListener('dragover', handleDragOver);
    fileUpload.addEventListener('dragleave', handleDragLeave);
    fileUpload.addEventListener('drop', handleDrop);
    fileInput.addEventListener('change', handleFileSelect);
}

function handleDragOver(e) {
    e.preventDefault();
    e.currentTarget.classList.add('dragover');
}

function handleDragLeave(e) {
    e.preventDefault();
    e.currentTarget.classList.remove('dragover');
}

function handleDrop(e) {
    e.preventDefault();
    e.currentTarget.classList.remove('dragover');
    const files = e.dataTransfer.files;
    if (files.length > 0) {
        processFile(files[0]);
    }
}

function handleFileSelect(e) {
    const files = e.target.files;
    if (files.length > 0) {
        processFile(files[0]);
    }
}

function processFile(file) {
    if (file.type !== 'text/csv' && !file.name.endsWith('.csv')) {
        showError('Please select a CSV file.');
        return;
    }

    Papa.parse(file, {
        complete: function(results) {
            if (results.errors.length > 0) {
                showError('Error parsing CSV: ' + results.errors[0].message);
                return;
            }

            // Extract numeric data from first column
            const data = results.data
                .map(row => parseFloat(row[0]))
                .filter(val => !isNaN(val));

            if (data.length === 0) {
                showError('No valid numeric data found in CSV file.');
                return;
            }

            currentData = data;
            document.getElementById('fileInfo').innerHTML = 
                `Loaded ${data.length} data points from ${file.name}`;
            
            showSuccess(`Successfully loaded ${data.length} data points.`);
            
            // Show raw data visualization
            showRawDataVisualization();
        },
        error: function(error) {
            showError('Error reading file: ' + error.message);
        }
    });
}

// Initialize WebAssembly module
function initializeWebAssembly() {
    // This would normally load the compiled WebAssembly module
    // For now, we'll simulate the functionality with JavaScript
    console.log('WebAssembly module would be loaded here');
    
    // Simulate module loading
    setTimeout(() => {
        Module = {
            _set_data: function(ptr, length) {
                console.log('Setting data:', length, 'points');
            },
            _set_parameters: function(iterations, burnin, thin) {
                console.log('Setting parameters:', iterations, burnin, thin);
            },
            _run_mcmc_analysis: function() {
                console.log('Running MCMC analysis...');
                return 0;
            },
            _get_results_size: function() {
                return 1000;
            },
            _get_component_counts: function() {
                return 0;
            }
        };
    }, 1000);
}

// Load demo data (simulated mixture of 3 Gaussians)
function loadDemoData() {
    const n = 1000;
    const data = [];
    
    // Generate data from mixture of 3 Gaussians
    // Component 1: N(2, 1) with weight 0.3
    for (let i = 0; i < n * 0.3; i++) {
        data.push(2 + (Math.random() - 0.5) * 4); // Approximate N(2,1)
    }
    
    // Component 2: N(6, 0.5) with weight 0.4
    for (let i = 0; i < n * 0.4; i++) {
        data.push(6 + (Math.random() - 0.5) * 2); // Approximate N(6,0.5)
    }
    
    // Component 3: N(12, 1) with weight 0.3
    for (let i = 0; i < n * 0.3; i++) {
        data.push(12 + (Math.random() - 0.5) * 4); // Approximate N(12,1)
    }
    
    // Shuffle the data
    for (let i = data.length - 1; i > 0; i--) {
        const j = Math.floor(Math.random() * (i + 1));
        [data[i], data[j]] = [data[j], data[i]];
    }
    
    currentData = data;
    document.getElementById('fileInfo').innerHTML = 
        `Loaded ${data.length} demo data points (3-component mixture)`;
    
    showSuccess(`Demo data loaded with ${data.length} points from a 3-component mixture.`);
    
    // Show raw data visualization
    showRawDataVisualization();
}

// Run RJMCMC analysis with dynamic visualization
function runAnalysis() {
    if (!currentData || currentData.length === 0) {
        showError('Please load data first.');
        return;
    }

    if (isRunning) {
        showError('Analysis is already running.');
        return;
    }

    // Get parameters
    const iterations = parseInt(document.getElementById('iterations').value);
    const burnin = parseInt(document.getElementById('burnin').value);
    const thin = parseInt(document.getElementById('thin').value);
    const maxComponents = parseInt(document.getElementById('maxComponents').value);

    totalIterations = iterations;
    currentIteration = 0;
    isRunning = true;

    // Show loading indicator
    document.getElementById('loadingIndicator').style.display = 'block';
    document.getElementById('resultsContainer').style.display = 'block';
    document.getElementById('welcomeMessage').style.display = 'none';
    document.getElementById('runAnalysis').disabled = true;
    
    // Show raw data visualization
    showRawDataVisualization();

    // Initialize MCMC state
    mcmcState = {
        nComponents: 2,
        components: [
            { p: 0.5, mu: 1.0, sigma: 1.0 },
            { p: 0.5, mu: 4.0, sigma: 1.0 }
        ],
        logPosterior: 0,
        trace: {
            components: [],
            logPosterior: [],
            iterations: []
        }
    };

    // Initialize charts
    initializeTraceCharts();
    
    // Start dynamic MCMC
    runDynamicMCMC(iterations, burnin, thin, maxComponents);
}

function runDynamicMCMC(iterations, burnin, thin, maxComponents) {
    const updateInterval = 50; // Update every 50 iterations
    let iteration = 0;
    
    const mcmcStep = () => {
        if (iteration >= iterations || !isRunning) {
            // Analysis complete
            isRunning = false;
            document.getElementById('runAnalysis').disabled = false;
            document.getElementById('loadingIndicator').style.display = 'none';
            updateFinalResults();
            return;
        }

        // Perform MCMC step
        performMCMCStep(maxComponents);
        
        // Update progress
        currentIteration = iteration;
        const progress = (iteration / iterations) * 100;
        document.getElementById('progressFill').style.width = progress + '%';
        document.getElementById('iterationInfo').textContent = 
            `Iteration: ${iteration} / ${iterations}`;

        // Update trace plots every few iterations
        if (iteration % 10 === 0) {
            updateTraceCharts();
        }

        // Update main charts less frequently
        if (iteration % 100 === 0) {
            updateMainCharts();
        }

        iteration++;
        
        // Continue with next step
        setTimeout(mcmcStep, 10); // Small delay for visualization
    };

    mcmcStep();
}

function performMCMCStep(maxComponents) {
    // Simulate RJMCMC step
    const moveType = Math.random();
    
    if (mcmcState.nComponents === 1) {
        if (moveType < 0.5) {
            birthMove();
        } else {
            jumpMove();
        }
    } else if (mcmcState.nComponents >= maxComponents) {
        if (moveType < 0.5) {
            deathMove();
        } else {
            jumpMove();
        }
    } else {
        if (moveType < 0.33) {
            deathMove();
        } else if (moveType < 0.67) {
            jumpMove();
        } else {
            birthMove();
        }
    }
    
    // Calculate log posterior (simplified)
    mcmcState.logPosterior = calculateLogPosterior();
    
    // Store trace
    mcmcState.trace.components.push(mcmcState.nComponents);
    mcmcState.trace.logPosterior.push(mcmcState.logPosterior);
    mcmcState.trace.iterations.push(currentIteration);
}

function birthMove() {
    if (mcmcState.nComponents >= 8) return;
    
    // Add new component
    const newComponent = {
        p: Math.random() * 0.3,
        mu: Math.random() * 20,
        sigma: 0.3 + Math.random() * 2.7
    };
    
    // Adjust existing proportions
    const totalP = mcmcState.components.reduce((sum, comp) => sum + comp.p, 0);
    const scaleFactor = (totalP - newComponent.p) / totalP;
    
    mcmcState.components.forEach(comp => {
        comp.p *= scaleFactor;
    });
    
    mcmcState.components.push(newComponent);
    mcmcState.nComponents++;
}

function deathMove() {
    if (mcmcState.nComponents <= 1) return;
    
    // Remove random component
    const idx = Math.floor(Math.random() * mcmcState.nComponents);
    const removedP = mcmcState.components[idx].p;
    
    mcmcState.components.splice(idx, 1);
    mcmcState.nComponents--;
    
    // Renormalize proportions
    const scaleFactor = 1 / (1 - removedP);
    mcmcState.components.forEach(comp => {
        comp.p *= scaleFactor;
    });
}

function jumpMove() {
    // Update random component
    const idx = Math.floor(Math.random() * mcmcState.nComponents);
    const comp = mcmcState.components[idx];
    
    // Small random changes
    comp.p = Math.max(0.01, Math.min(0.99, comp.p + (Math.random() - 0.5) * 0.1));
    comp.mu = Math.max(0, Math.min(20, comp.mu + (Math.random() - 0.5) * 0.5));
    comp.sigma = Math.max(0.3, Math.min(3, comp.sigma + (Math.random() - 0.5) * 0.2));
    
    // Renormalize proportions
    const totalP = mcmcState.components.reduce((sum, c) => sum + c.p, 0);
    mcmcState.components.forEach(c => {
        c.p /= totalP;
    });
}

function calculateLogPosterior() {
    // Simplified log posterior calculation
    let logLik = 0;
    
    // Sample a few data points for likelihood
    const sampleSize = Math.min(100, currentData.length);
    for (let i = 0; i < sampleSize; i++) {
        const idx = Math.floor(Math.random() * currentData.length);
        const x = currentData[idx];
        
        let mixtureProb = 0;
        mcmcState.components.forEach(comp => {
            mixtureProb += comp.p * normalDensity(x, comp.mu, comp.sigma);
        });
        
        if (mixtureProb > 0) {
            logLik += Math.log(mixtureProb);
        }
    }
    
    // Prior (simplified)
    const logPrior = -Math.log(8) - mcmcState.nComponents * 2; // Penalty for more components
    
    return logLik + logPrior;
}

function normalDensity(x, mu, sigma) {
    return (1 / (sigma * Math.sqrt(2 * Math.PI))) * 
           Math.exp(-0.5 * Math.pow((x - mu) / sigma, 2));
}

function initializeTraceCharts() {
    // Component trace chart
    const componentCtx = document.getElementById('componentTraceChart').getContext('2d');
    componentTraceChart = new Chart(componentCtx, {
        type: 'line',
        data: {
            labels: [],
            datasets: [{
                label: 'Components',
                data: [],
                borderColor: '#000000',
                backgroundColor: 'rgba(0,0,0,0.1)',
                borderWidth: 1,
                fill: false,
                tension: 0
            }]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            scales: {
                y: {
                    beginAtZero: true,
                    ticks: {
                        stepSize: 1
                    }
                },
                x: {
                    display: false
                }
            },
            plugins: {
                legend: {
                    display: false
                }
            }
        }
    });

    // Log posterior trace chart
    const logPostCtx = document.getElementById('logPosteriorChart').getContext('2d');
    logPosteriorChart = new Chart(logPostCtx, {
        type: 'line',
        data: {
            labels: [],
            datasets: [{
                label: 'Log Posterior',
                data: [],
                borderColor: '#000000',
                backgroundColor: 'rgba(0,0,0,0.1)',
                borderWidth: 1,
                fill: false,
                tension: 0
            }]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            scales: {
                y: {
                    beginAtZero: false,
                    display: true
                },
                x: {
                    display: false
                }
            },
            plugins: {
                legend: {
                    display: false
                }
            }
        }
    });
}

function updateTraceCharts() {
    if (!mcmcState || !componentTraceChart || !logPosteriorChart) return;
    
    const trace = mcmcState.trace;
    const isRunning = mcmcState.isRunning || false;
    const showBurnin = mcmcState.showBurnin !== false; // Default to true
    
    // Helper function to get data range (same as in rjmc_simple.js)
    function getDataRange(data, isRunning = false, showBurnin = true) {
        if (!data || data.length === 0) return { data: [], labels: [] };
        
        if (isRunning) {
            // During run: show only last 1000 points for performance
            const maxPoints = 1000;
            const startIndex = Math.max(0, data.length - maxPoints);
            return {
                data: data.slice(startIndex),
                labels: data.slice(startIndex).map((_, i) => startIndex + i)
            };
        } else {
            // After run: show full trajectory
            if (showBurnin) {
                return {
                    data: data,
                    labels: data.map((_, i) => i)
                };
            } else {
                // Remove burn-in samples
                const burnin = Math.min(100, data.length);
                return {
                    data: data.slice(burnin),
                    labels: data.slice(burnin).map((_, i) => burnin + i)
                };
            }
        }
    }
    
    // Update component trace with smart data range
    const { data: componentData, labels: componentLabels } = getDataRange(
        trace.components, isRunning, showBurnin
    );
    
    componentTraceChart.data.labels = componentLabels;
    componentTraceChart.data.datasets[0].data = componentData;
    componentTraceChart.update('none');
    
    // Update log posterior trace with smart data range
    const { data: logPostData, labels: logPostLabels } = getDataRange(
        trace.logPosterior, isRunning, showBurnin
    );
    
    logPosteriorChart.data.labels = logPostLabels;
    logPosteriorChart.data.datasets[0].data = logPostData;
    logPosteriorChart.update('none');
}

function updateMainCharts() {
    if (!mcmcState) return;
    
    // Update component distribution
    updateComponentChart();
    
    // Update fit chart
    updateFitChart();
}

function updateComponentChart() {
    if (!mcmcState) return;
    
    const ctx = document.getElementById('componentChart').getContext('2d');
    
    if (componentChart) {
        componentChart.destroy();
    }
    
    // Calculate component counts from trace
    const componentCounts = {};
    mcmcState.trace.components.forEach(n => {
        componentCounts[n] = (componentCounts[n] || 0) + 1;
    });
    
    const labels = Object.keys(componentCounts).sort((a, b) => a - b);
    const data = labels.map(label => componentCounts[label] / mcmcState.trace.components.length);
    
    componentChart = new Chart(ctx, {
        type: 'bar',
        data: {
            labels: labels,
            datasets: [{
                label: 'Posterior Probability',
                data: data,
                backgroundColor: '#000000',
                borderColor: '#000000',
                borderWidth: 1
            }]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            scales: {
                y: {
                    beginAtZero: true,
                    max: Math.max(...data) * 1.1, // Dynamic max with 10% padding
                    title: {
                        display: true,
                        text: 'Probability'
                    }
                },
                x: {
                    title: {
                        display: true,
                        text: 'Number of Components'
                    }
                }
            },
            plugins: {
                legend: {
                    display: false
                }
            }
        }
    });
}

function updateFitChart() {
    if (!mcmcState || !currentData) return;
    
    const ctx = document.getElementById('fitChart').getContext('2d');
    
    if (fitChart) {
        fitChart.destroy();
    }
    
    // Create histogram of data
    const histogram = createHistogram(currentData, 30);
    
    // Create fitted mixture model
    const fittedCurve = createFittedCurve(histogram.x, mcmcState.components);
    
    fitChart = new Chart(ctx, {
        type: 'line',
        data: {
            labels: histogram.x,
            datasets: [
                {
                    label: 'Observed Data',
                    data: histogram.y,
                    type: 'bar',
                    backgroundColor: 'rgba(0,0,0,0.3)',
                    borderColor: '#000000',
                    borderWidth: 1
                },
                {
                    label: 'Fitted Mixture',
                    data: fittedCurve,
                    type: 'line',
                    borderColor: '#000000',
                    backgroundColor: 'rgba(0,0,0,0.1)',
                    borderWidth: 2,
                    fill: true
                }
            ]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            scales: {
                y: {
                    beginAtZero: true,
                    title: {
                        display: true,
                        text: 'Density'
                    }
                },
                x: {
                    title: {
                        display: true,
                        text: 'Value'
                    }
                }
            }
        }
    });
}

function createHistogram(data, bins) {
    const min = Math.min(...data);
    const max = Math.max(...data);
    const binWidth = (max - min) / bins;
    
    const histogram = Array(bins).fill(0);
    const x = [];
    
    for (let i = 0; i < bins; i++) {
        x.push(min + i * binWidth + binWidth / 2);
    }
    
    data.forEach(value => {
        const binIndex = Math.min(Math.floor((value - min) / binWidth), bins - 1);
        histogram[binIndex]++;
    });
    
    // Normalize
    const total = data.length;
    const y = histogram.map(count => count / (total * binWidth));
    
    return { x, y };
}

function createFittedCurve(x, components) {
    return x.map(xVal => {
        let density = 0;
        components.forEach(comp => {
            density += comp.p * normalDensity(xVal, comp.mu, comp.sigma);
        });
        return density;
    });
}

function showRawDataVisualization() {
    if (!currentData || currentData.length === 0) return;
    
    const ctx = document.getElementById('rawDataChart').getContext('2d');
    
    if (rawDataChart) {
        rawDataChart.destroy();
    }
    
    // Create histogram of raw data
    const histogram = createHistogram(currentData, 40);
    
    // Calculate basic statistics
    const mean = currentData.reduce((sum, val) => sum + val, 0) / currentData.length;
    const variance = currentData.reduce((sum, val) => sum + Math.pow(val - mean, 2), 0) / currentData.length;
    const stdDev = Math.sqrt(variance);
    const min = Math.min(...currentData);
    const max = Math.max(...currentData);
    
    rawDataChart = new Chart(ctx, {
        type: 'bar',
        data: {
            labels: histogram.x,
            datasets: [{
                label: 'Data Frequency',
                data: histogram.y,
                backgroundColor: 'rgba(0,0,0,0.3)',
                borderColor: '#000000',
                borderWidth: 1
            }]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            scales: {
                y: {
                    beginAtZero: true,
                    title: {
                        display: true,
                        text: 'Frequency'
                    }
                },
                x: {
                    title: {
                        display: true,
                        text: 'Value'
                    }
                }
            },
            plugins: {
                legend: {
                    display: false
                },
                tooltip: {
                    callbacks: {
                        title: function(context) {
                            const value = context[0].label;
                            return `Value: ${parseFloat(value).toFixed(2)}`;
                        },
                        label: function(context) {
                            const count = context.parsed.y;
                            const total = currentData.length;
                            const percentage = ((count / total) * 100).toFixed(1);
                            return `Count: ${count} (${percentage}%)`;
                        }
                    }
                }
            }
        }
    });
    
    // Add data statistics to the chart container
    const chartContainer = document.querySelector('#rawDataChart').closest('.chart-container');
    let statsDiv = chartContainer.querySelector('.data-stats');
    
    if (!statsDiv) {
        statsDiv = document.createElement('div');
        statsDiv.className = 'data-stats';
        statsDiv.style.cssText = `
            margin-top: 15px;
            padding: 15px;
            border: 1px solid #000000;
            font-size: 0.8em;
            background: #f8f8f8;
        `;
        chartContainer.appendChild(statsDiv);
    }
    
    statsDiv.innerHTML = `
        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(120px, 1fr)); gap: 10px;">
            <div><strong>Count:</strong> ${currentData.length}</div>
            <div><strong>Mean:</strong> ${mean.toFixed(3)}</div>
            <div><strong>Std Dev:</strong> ${stdDev.toFixed(3)}</div>
            <div><strong>Min:</strong> ${min.toFixed(3)}</div>
            <div><strong>Max:</strong> ${max.toFixed(3)}</div>
            <div><strong>Range:</strong> ${(max - min).toFixed(3)}</div>
        </div>
    `;
}

function updateFinalResults() {
    if (!mcmcState) return;
    
    // Calculate final statistics
    const componentCounts = {};
    mcmcState.trace.components.forEach(n => {
        componentCounts[n] = (componentCounts[n] || 0) + 1;
    });
    
    const mostLikelyComponents = Object.keys(componentCounts).reduce((a, b) => 
        componentCounts[a] > componentCounts[b] ? a : b
    );
    
    const totalSamples = mcmcState.trace.components.length;
    const probability3Components = (componentCounts[3] || 0) / totalSamples;
    
    // Update summary stats
    document.getElementById('summaryStats').innerHTML = `
        <div class="stat-item">
            <div class="stat-value">${mostLikelyComponents}</div>
            <div class="stat-label">Most Likely Components</div>
        </div>
        <div class="stat-item">
            <div class="stat-value">${totalSamples}</div>
            <div class="stat-label">Posterior Samples</div>
        </div>
        <div class="stat-item">
            <div class="stat-value">${currentData.length}</div>
            <div class="stat-label">Data Points</div>
        </div>
        <div class="stat-item">
            <div class="stat-value">${probability3Components.toFixed(3)}</div>
            <div class="stat-label">P(K=3)</div>
        </div>
    `;
}

// Utility functions
function showError(message) {
    const errorDiv = document.createElement('div');
    errorDiv.className = 'error';
    errorDiv.textContent = message;
    
    const container = document.querySelector('.visualization-panel');
    container.insertBefore(errorDiv, container.firstChild);
    
    setTimeout(() => errorDiv.remove(), 5000);
}

function showSuccess(message) {
    const successDiv = document.createElement('div');
    successDiv.className = 'success';
    successDiv.textContent = message;
    
    const container = document.querySelector('.visualization-panel');
    container.insertBefore(successDiv, container.firstChild);
    
    setTimeout(() => successDiv.remove(), 3000);
}