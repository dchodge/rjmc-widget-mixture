// R WebAssembly Application for RJMCMC Widget
// This replaces the JavaScript simulation with actual R WebAssembly calls

let rWebModule = null;
let currentData = null;
let analysisResults = null;
let isAnalysisRunning = false;

// Chart instances for dynamic updates
let rawDataChart = null;
let componentTraceChart = null;
let logPosteriorChart = null;
let componentDistChart = null;
let mixtureFitChart = null;

// Initialize the R WebAssembly module
async function initializeRWebModule() {
    try {
        console.log('[R-WASM] Initializing R WebAssembly module...');
        
        // Load R WebAssembly module
        rWebModule = await webr.install();
        console.log('[R-WASM] R WebAssembly module loaded successfully');
        
        // Load required packages
        await rWebModule.webR.evalR(`
            if (!require("webr", quietly = TRUE)) {
                install.packages("webr", repos = "https://webr.r-wasm.org/latest")
            }
        `);
        
        // Load the rjmc package
        await rWebModule.webR.evalR(`
            if (!require("rjmc", quietly = TRUE)) {
                devtools::load_all("/Users/davidhodgson/Dropbox/Mac (3)/Documents/research/software/mcmc/rjmc")
            }
        `);
        
        // Load the RJMCMC server functions
        await rWebModule.webR.evalR(`
            source('r_wasm_server.R')
        `);
        
        console.log('[R-WASM] RJMCMC server loaded successfully');
        return true;
        
    } catch (error) {
        console.error('[R-WASM] Failed to initialize R WebAssembly module:', error);
        return false;
    }
}

// Generate sample data using R
async function generateSampleData() {
    try {
        console.log('[R-WASM] Generating sample data...');
        
        const result = await rWebModule.webR.evalR(`
            generate_sample_data(1000)
        `);
        
        const data = await result.toJs();
        currentData = data;
        
        console.log('[R-WASM] Sample data generated:', data.obs.length, 'points');
        return data;
        
    } catch (error) {
        console.error('[R-WASM] Failed to generate sample data:', error);
        throw error;
    }
}

// Run RJMCMC analysis with streaming updates
async function runRJMCMCAnalysis(settings = {}) {
    if (isAnalysisRunning) {
        console.log('[R-WASM] Analysis already running, skipping...');
        return;
    }
    
    if (!currentData || !currentData.obs || currentData.obs.length === 0) {
        throw new Error('No data available for analysis');
    }
    
    isAnalysisRunning = true;
    
    try {
        console.log('[R-WASM] Starting RJMCMC analysis...');
        
        // Default settings
        const defaultSettings = {
            numberCores: 1,
            numberChainRuns: 4,
            iterations: 20000,
            burninPosterior: 10000,
            thin: 5,
            runParallel: false
        };
        
        const finalSettings = { ...defaultSettings, ...settings };
        
        // Show initial progress
        updateProgress(0, finalSettings.iterations, 'Initializing RJMCMC...');
        
        // Run RJMCMC
        const result = await rWebModule.webR.evalR(`
            run_rjmc_streaming(
                list(obs = ${JSON.stringify(currentData.obs)}),
                ${JSON.stringify(finalSettings)}
            )
        `);
        
        const analysisResult = await result.toJs();
        
        if (analysisResult.success) {
            analysisResults = analysisResult;
            console.log('[R-WASM] RJMCMC analysis completed successfully');
            
            // Get posterior summary
            const summaryResult = await rWebModule.webR.evalR(`
                get_posterior_summary(outputs)
            `);
            
            const summary = await summaryResult.toJs();
            
            if (summary.success) {
                console.log('[R-WASM] Most probable K:', summary.most_probable_K);
                updateAnalysisResults(summary);
            }
            
            updateProgress(finalSettings.iterations, finalSettings.iterations, 'Analysis complete!');
            
        } else {
            throw new Error(analysisResult.error || 'RJMCMC analysis failed');
        }
        
    } catch (error) {
        console.error('[R-WASM] RJMCMC analysis failed:', error);
        updateProgress(0, 0, 'Analysis failed: ' + error.message);
        throw error;
        
    } finally {
        isAnalysisRunning = false;
    }
}

// Update progress display
function updateProgress(current, total, message) {
    const progressElement = document.getElementById('mcmc-progress');
    const progressBar = document.getElementById('progress-bar');
    
    if (progressElement && progressBar) {
        const percentage = total > 0 ? (current / total * 100).toFixed(1) : 0;
        progressElement.textContent = `Step ${current} / ${total} (${percentage}%) - ${message}`;
        progressBar.style.width = percentage + '%';
    }
    
    console.log(`[PROGRESS] ${current}/${total} (${percentage}%) - ${message}`);
}

// Update analysis results and visualizations
function updateAnalysisResults(summary) {
    try {
        // Update component count display
        const kElement = document.getElementById('current-k');
        if (kElement) {
            kElement.textContent = summary.most_probable_K;
        }
        
        // Update K probabilities
        updateKProbabilities(summary.K_probabilities);
        
        // Update visualizations
        updateAllCharts(summary);
        
        console.log('[R-WASM] Analysis results updated successfully');
        
    } catch (error) {
        console.error('[R-WASM] Failed to update analysis results:', error);
    }
}

// Update K probabilities display
function updateKProbabilities(kProbs) {
    const kProbsElement = document.getElementById('k-probabilities');
    if (kProbsElement) {
        let html = '<h4>Component Probabilities:</h4>';
        for (const [k, prob] of Object.entries(kProbs)) {
            const percentage = (prob * 100).toFixed(1);
            html += `<div>K=${k}: ${percentage}%</div>`;
        }
        kProbsElement.innerHTML = html;
    }
}

// Update all charts with new data
function updateAllCharts(summary) {
    try {
        // Update component trace plot
        updateComponentTracePlot();
        
        // Update log posterior plot
        updateLogPosteriorPlot();
        
        // Update component distribution
        updateComponentDistribution(summary);
        
        // Update mixture fit
        updateMixtureFit(summary);
        
    } catch (error) {
        console.error('[R-WASM] Failed to update charts:', error);
    }
}

// Update component trace plot
function updateComponentTracePlot() {
    if (!analysisResults || !analysisResults.outputs) return;
    
    try {
        const ctx = document.getElementById('componentTraceChart').getContext('2d');
        
        if (componentTraceChart) {
            componentTraceChart.destroy();
        }
        
        // Extract K values from outputs
        const kValues = [];
        const iterations = [];
        
        // This would need to be implemented based on the actual RJMCMC output structure
        // For now, create a placeholder
        for (let i = 0; i < 100; i++) {
            iterations.push(i);
            kValues.push(2 + Math.random() * 2); // Placeholder
        }
        
        componentTraceChart = new Chart(ctx, {
            type: 'line',
            data: {
                labels: iterations,
                datasets: [{
                    label: 'Number of Components (K)',
                    data: kValues,
                    borderColor: '#000000',
                    backgroundColor: 'rgba(0,0,0,0.1)',
                    borderWidth: 2,
                    fill: false
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
                            text: 'Number of Components'
                        }
                    },
                    x: {
                        title: {
                            display: true,
                            text: 'Iteration'
                        }
                    }
                }
            }
        });
        
    } catch (error) {
        console.error('[R-WASM] Failed to update component trace plot:', error);
    }
}

// Update log posterior plot
function updateLogPosteriorPlot() {
    if (!analysisResults || !analysisResults.outputs) return;
    
    try {
        const ctx = document.getElementById('logPosteriorChart').getContext('2d');
        
        if (logPosteriorChart) {
            logPosteriorChart.destroy();
        }
        
        // Extract log posterior values
        const logPostValues = [];
        const iterations = [];
        
        // Placeholder data - would need actual implementation
        for (let i = 0; i < 100; i++) {
            iterations.push(i);
            logPostValues.push(-1000 + Math.random() * 100); // Placeholder
        }
        
        logPosteriorChart = new Chart(ctx, {
            type: 'line',
            data: {
                labels: iterations,
                datasets: [{
                    label: 'Log Posterior',
                    data: logPostValues,
                    borderColor: '#000000',
                    backgroundColor: 'rgba(0,0,0,0.1)',
                    borderWidth: 2,
                    fill: false
                }]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                scales: {
                    y: {
                        title: {
                            display: true,
                            text: 'Log Posterior'
                        }
                    },
                    x: {
                        title: {
                            display: true,
                            text: 'Iteration'
                        }
                    }
                }
            }
        });
        
    } catch (error) {
        console.error('[R-WASM] Failed to update log posterior plot:', error);
    }
}

// Update component distribution
function updateComponentDistribution(summary) {
    if (!summary || !summary.posterior_samples) return;
    
    try {
        const ctx = document.getElementById('componentDistChart').getContext('2d');
        
        if (componentDistChart) {
            componentDistChart.destroy();
        }
        
        // Create histogram of component counts
        const kCounts = {};
        for (const [k, prob] of Object.entries(summary.K_probabilities)) {
            kCounts[k] = prob * 100; // Convert to percentage
        }
        
        componentDistChart = new Chart(ctx, {
            type: 'bar',
            data: {
                labels: Object.keys(kCounts),
                datasets: [{
                    label: 'Probability (%)',
                    data: Object.values(kCounts),
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
                            text: 'Probability (%)'
                        }
                    },
                    x: {
                        title: {
                            display: true,
                            text: 'Number of Components (K)'
                        }
                    }
                }
            }
        });
        
    } catch (error) {
        console.error('[R-WASM] Failed to update component distribution:', error);
    }
}

// Update mixture fit
function updateMixtureFit(summary) {
    if (!currentData || !summary) return;
    
    try {
        const ctx = document.getElementById('mixtureFitChart').getContext('2d');
        
        if (mixtureFitChart) {
            mixtureFitChart.destroy();
        }
        
        // Create histogram of raw data
        const histogram = createHistogram(currentData.obs, 40);
        
        // Create fitted mixture density (placeholder)
        const fittedDensity = histogram.x.map(x => {
            // This would be calculated from the actual RJMCMC results
            return Math.exp(-0.5 * Math.pow((x - 6) / 2, 2)) / (2 * Math.sqrt(2 * Math.PI));
        });
        
        mixtureFitChart = new Chart(ctx, {
            type: 'line',
            data: {
                labels: histogram.x,
                datasets: [
                    {
                        label: 'Data',
                        data: histogram.y,
                        type: 'bar',
                        backgroundColor: 'rgba(0,0,0,0.3)',
                        borderColor: '#000000',
                        borderWidth: 1
                    },
                    {
                        label: 'Fitted Mixture',
                        data: fittedDensity,
                        type: 'line',
                        borderColor: '#000000',
                        backgroundColor: 'transparent',
                        borderWidth: 2,
                        fill: false
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
        
    } catch (error) {
        console.error('[R-WASM] Failed to update mixture fit:', error);
    }
}

// Utility function to create histogram
function createHistogram(data, bins) {
    const min = Math.min(...data);
    const max = Math.max(...data);
    const binWidth = (max - min) / bins;
    
    const histogram = {
        x: [],
        y: []
    };
    
    for (let i = 0; i < bins; i++) {
        const binStart = min + i * binWidth;
        const binEnd = min + (i + 1) * binWidth;
        const binCenter = (binStart + binEnd) / 2;
        
        const count = data.filter(d => d >= binStart && d < binEnd).length;
        
        histogram.x.push(binCenter.toFixed(2));
        histogram.y.push(count);
    }
    
    return histogram;
}

// Show raw data visualization
function showRawDataVisualization() {
    if (!currentData || currentData.obs.length === 0) return;
    
    const ctx = document.getElementById('rawDataChart').getContext('2d');
    
    if (rawDataChart) {
        rawDataChart.destroy();
    }
    
    // Create histogram of raw data
    const histogram = createHistogram(currentData.obs, 40);
    
    // Calculate basic statistics
    const mean = currentData.obs.reduce((sum, val) => sum + val, 0) / currentData.obs.length;
    const variance = currentData.obs.reduce((sum, val) => sum + Math.pow(val - mean, 2), 0) / currentData.obs.length;
    const stdDev = Math.sqrt(variance);
    const min = Math.min(...currentData.obs);
    const max = Math.max(...currentData.obs);
    
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
            <div><strong>Count:</strong> ${currentData.obs.length}</div>
            <div><strong>Mean:</strong> ${mean.toFixed(3)}</div>
            <div><strong>Std Dev:</strong> ${stdDev.toFixed(3)}</div>
            <div><strong>Min:</strong> ${min.toFixed(3)}</div>
            <div><strong>Max:</strong> ${max.toFixed(3)}</div>
            <div><strong>Range:</strong> ${(max - min).toFixed(3)}</div>
        </div>
    `;
}

// Initialize the application
async function initRWebApp() {
    try {
        console.log('[R-WASM] Initializing R WebAssembly application...');
        
        // Initialize R WebAssembly module
        const success = await initializeRWebModule();
        if (!success) {
            throw new Error('Failed to initialize R WebAssembly module');
        }
        
        // Generate initial sample data
        await generateSampleData();
        
        // Show raw data visualization
        showRawDataVisualization();
        
        // Enable analysis button
        document.getElementById('mcmc-btn').disabled = false;
        document.getElementById('mcmc-btn').textContent = 'Start RJMCMC Analysis';
        
        console.log('[R-WASM] Application initialized successfully');
        
    } catch (error) {
        console.error('[R-WASM] Failed to initialize application:', error);
        document.getElementById('error-message').style.display = 'block';
        document.getElementById('error-text').textContent = `Initialization Error: ${error.message}`;
    }
}

// Start RJMCMC analysis
async function startRJMCMCAnalysis() {
    try {
        console.log('[R-WASM] Starting RJMCMC analysis...');
        
        // Get settings from UI
        const iterations = parseInt(document.getElementById('mcmc_steps').value);
        const burnin = parseInt(document.getElementById('burnin').value);
        
        const settings = {
            iterations: iterations,
            burninPosterior: burnin,
            numberChainRuns: 4,
            numberCores: 1,
            thin: 5,
            runParallel: false
        };
        
        // Disable button during analysis
        document.getElementById('mcmc-btn').disabled = true;
        document.getElementById('mcmc-btn').textContent = 'Running RJMCMC...';
        
        // Run analysis
        await runRJMCMCAnalysis(settings);
        
        // Re-enable button
        document.getElementById('mcmc-btn').disabled = false;
        document.getElementById('mcmc-btn').textContent = 'Start RJMCMC Analysis';
        
    } catch (error) {
        console.error('[R-WASM] RJMCMC analysis failed:', error);
        
        // Re-enable button
        document.getElementById('mcmc-btn').disabled = false;
        document.getElementById('mcmc-btn').textContent = 'Start RJMCMC Analysis';
        
        // Show error
        document.getElementById('error-message').style.display = 'block';
        document.getElementById('error-text').textContent = `Analysis Error: ${error.message}`;
    }
}

// Make functions globally available
window.initRWebApp = initRWebApp;
window.startRJMCMCAnalysis = startRJMCMCAnalysis;
window.showRawDataVisualization = showRawDataVisualization;

// Initialize when DOM is loaded
document.addEventListener('DOMContentLoaded', initRWebApp);
