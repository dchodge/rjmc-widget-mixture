// Simplified RJMCMC Implementation for Web Widget
// This implements the core RJMCMC algorithm from the vignette in JavaScript
// for immediate functionality while R WebAssembly is being set up

let currentData = null;
let analysisResults = null;
let isAnalysisRunning = false;

// Prior configuration for K regularization
const PRIOR_CONFIG = {
    // Poisson prior parameters
    poissonLambda: 2.0,        // Mean of Poisson prior (favors K around 2)
    
    // Additional regularization
    regularizationStrength: 0.1, // Weak quadratic penalty on K
    
    // Prior bounds
    minK: 1,
    maxK: 8
};

// Chart instances for dynamic updates
// Raw data chart now uses Plotly (no variable needed)
// Component trace chart now uses Plotly (no variable needed)
let logPosteriorChart = null;
let componentDistChart = null;
let mixtureFitChart = null;
let meanTraceChart = null;
let sigmaTraceChart = null;

// RJMCMC state
let rjmcState = {
    chains: [], // Array of 4 chains
    currentK: 2,
    currentParams: { sigma: 0 },
    currentJump: null,
    logPosterior: -Infinity,
    iterations: 0,
    acceptedMoves: 0,
    kHistory: [], // Combined history from all chains
    logPostHistory: [], // Combined history from all chains
    kProbabilities: {},
    numChains: 4,
    topKFits: {}, // Store best fits for each K
    currentTab: 'top1',
    kPosteriorSamples: {}, // Store all posterior samples for each K
    kPosteriorMeans: {}, // Store posterior means for each K
    isAnalysisRunning: false, // Track if analysis is running
    parameterTraces: {
        mu: [], // Store mean traces for each chain
        sigma: [] // Store sigma traces for each chain
    },
    maxHistoryLength: 50000, // Increased limit for full chain display
    isRunning: false, // Track if MCMC is currently running
    showBurnin: true, // Toggle for showing burn-in samples
    dynamicParameterPlots: false, // Toggle for dynamic parameter plot updates
    updateFrequency: 10, // Update plots every N iterations (will be dynamic)
    cleanupFrequency: 500, // Clean up old data every N iterations (more frequent)
    maxDisplayPoints: 1000 // Maximum points to display during run
};

// Helper function to get data range based on MCMC state
function getDataRange(data, isRunning = false, showBurnin = true) {
    if (!data || data.length === 0) return { data: [], labels: [] };
    
    if (isRunning) {
        // During run: show only last N points for performance (adaptive)
        const maxPoints = rjmcState.maxDisplayPoints || 1000;
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

// Generate sample data like in the vignette
function generateSampleData() {
    const m = 5;
    const obs = [];
    
    // Generate data from 3-component mixture: K=3, μ=[2,6,12], σ=[1,0.5,1], π=[0.3,0.7,1]
    for (let i = 0; i < 30 * m; i++) {
        obs.push(2 + gaussianRandom() * 1); // Component 1: μ=2, σ=1
    }
    for (let i = 0; i < 70 * m; i++) {
        obs.push(6 + gaussianRandom() * 0.5); // Component 2: μ=6, σ=0.5
    }
    for (let i = 0; i < 100 * m; i++) {
        obs.push(12 + gaussianRandom() * 1); // Component 3: μ=12, σ=1
    }
    
    // Shuffle the data
    for (let i = obs.length - 1; i > 0; i--) {
        const j = Math.floor(Math.random() * (i + 1));
        [obs[i], obs[j]] = [obs[j], obs[i]];
    }
    
    return { obs: obs };
}

// Gaussian random number generator (Box-Muller transform)
function gaussianRandom() {
    if (gaussianRandom.spare !== undefined) {
        const temp = gaussianRandom.spare;
        delete gaussianRandom.spare;
        return temp;
    }
    
    const u1 = Math.random();
    const u2 = Math.random();
    const mag = Math.sqrt(-2 * Math.log(u1));
    gaussianRandom.spare = mag * Math.sin(2 * Math.PI * u2);
    return mag * Math.cos(2 * Math.PI * u2);
}

// Initialize RJMCMC state
function initializeRJMCMC() {
    // Initialize multiple chains
    rjmcState.chains = [];
    rjmcState.parameterTraces.mu = [];
    rjmcState.parameterTraces.sigma = [];
    
    for (let chainId = 0; chainId < rjmcState.numChains; chainId++) {
        const chain = {
            id: chainId,
            currentK: 2,
            currentParams: { sigma: gaussianRandom() },
            currentJump: initializeJumpMatrix(),
            logPosterior: -Infinity,
            kHistory: [],
            logPostHistory: [],
            acceptedMoves: 0,
            muHistory: [],
            sigmaHistory: []
        };
        chain.logPosterior = evaluateLogPosterior(chain.currentParams, chain.currentJump);
        rjmcState.chains.push(chain);
        
        // Initialize parameter trace arrays
        rjmcState.parameterTraces.mu.push([]);
        rjmcState.parameterTraces.sigma.push([]);
    }
    
    // Initialize combined state
    rjmcState.currentK = 2;
    rjmcState.currentParams = { sigma: gaussianRandom() };
    rjmcState.currentJump = initializeJumpMatrix();
    rjmcState.logPosterior = evaluateLogPosterior(rjmcState.currentParams, rjmcState.currentJump);
    rjmcState.iterations = 0;
    rjmcState.acceptedMoves = 0;
    rjmcState.kHistory = [];
    rjmcState.logPostHistory = [];
    rjmcState.kProbabilities = {};
    rjmcState.kPosteriorSamples = {};
    rjmcState.kPosteriorMeans = {};
    rjmcState.parameterTraces = {
        mu: [],
        sigma: []
    };
}

// Initialize jump matrix (mixture parameters)
function initializeJumpMatrix() {
    const p = [0.5, 0.5];
    const mu = [1, 4];
    const sigma = [1, 1];
    return { p: p, mu: mu, sigma: sigma };
}

// Log factorial function for Poisson prior
function logFactorial(n) {
    if (n <= 1) return 0;
    let result = 0;
    for (let i = 2; i <= n; i++) {
        result += Math.log(i);
    }
    return result;
}

// Evaluate log prior
function evaluateLogPrior(params, jump) {
    let logPrior = 0;
    
    // Prior for sigma
    logPrior += -0.5 * Math.pow(params.sigma, 2); // N(0,1)
    
    // Priors for mixture parameters
    for (let i = 0; i < jump.p.length; i++) {
        logPrior += Math.log(1); // Uniform(0,1) for p
        logPrior += Math.log(1/20); // Uniform(0,20) for mu
        logPrior += Math.log(1/2.7); // Uniform(0.3,3) for sigma
    }
    
    // Prior for K (number of components) - Poisson prior favoring smaller K
    const k = jump.p.length;
    const lambda = PRIOR_CONFIG.poissonLambda; // Poisson parameter (mean), favors K around 2
    logPrior += k * Math.log(lambda) - lambda - logFactorial(k);
    
    // Additional weak regularization term to further favor smaller K
    logPrior -= PRIOR_CONFIG.regularizationStrength * k * k; // Quadratic penalty
    
    return logPrior;
}

// Evaluate log likelihood
function evaluateLogLikelihood(params, jump) {
    let logLik = 0;
    const data = currentData.obs;
    
    for (let i = 0; i < data.length; i++) {
        let mixtureDensity = 0;
        for (let j = 0; j < jump.p.length; j++) {
            mixtureDensity += jump.p[j] * normalDensity(data[i], jump.mu[j], jump.sigma[j]);
        }
        logLik += Math.log(mixtureDensity);
    }
    
    return logLik;
}

// Normal density function
function normalDensity(x, mu, sigma) {
    return Math.exp(-0.5 * Math.pow((x - mu) / sigma, 2)) / (sigma * Math.sqrt(2 * Math.PI));
}

// Evaluate log posterior
function evaluateLogPosterior(params, jump) {
    const logPrior = evaluateLogPrior(params, jump);
    const logLik = evaluateLogLikelihood(params, jump);
    return logPrior + logLik;
}

// Sample birth proposal
function sampleBirthProposal(jump) {
    const pNew = Math.random();
    const pNewAdjusted = jump.p.map(p => p * (1 - pNew));
    pNewAdjusted.push(pNew);
    
    const muNew = [...jump.mu, Math.random() * 20];
    const sigmaNew = [...jump.sigma, 0.3 + Math.random() * 2.7];
    
    return { p: pNewAdjusted, mu: muNew, sigma: sigmaNew };
}

// Sample death proposal
function sampleDeathProposal(jump, indexToRemove) {
    const pNew = jump.p.filter((_, i) => i !== indexToRemove);
    const muNew = jump.mu.filter((_, i) => i !== indexToRemove);
    const sigmaNew = jump.sigma.filter((_, i) => i !== indexToRemove);
    
    // Renormalize probabilities
    const totalP = pNew.reduce((sum, p) => sum + p, 0);
    const pNormalized = pNew.map(p => p / totalP);
    
    return { p: pNormalized, mu: muNew, sigma: sigmaNew };
}

// Sample jump proposal (update existing component)
function sampleJumpProposal(jump, indexToUpdate) {
    const newJump = {
        p: [...jump.p],
        mu: [...jump.mu],
        sigma: [...jump.sigma]
    };
    
    // Update mixing proportion
    const pChange = gaussianRandom() * 0.01;
    newJump.p[indexToUpdate] = Math.max(0, Math.min(1, newJump.p[indexToUpdate] + pChange));
    
    // Renormalize all probabilities
    const totalP = newJump.p.reduce((sum, p) => sum + p, 0);
    newJump.p = newJump.p.map(p => p / totalP);
    
    // Update mean
    newJump.mu[indexToUpdate] += gaussianRandom() * 0.1;
    
    // Update standard deviation
    newJump.sigma[indexToUpdate] = Math.max(0.3, newJump.sigma[indexToUpdate] + gaussianRandom() * 0.1);
    
    return newJump;
}

// Sample proposal type
function sampleProposalType(currentK) {
    if (currentK === 2) {
        return Math.random() < 0.67 ? 'jump' : 'birth';
    } else if (currentK >= 8) {
        return Math.random() < 0.33 ? 'death' : 'jump';
    } else {
        const rand = Math.random();
        if (rand < 0.33) return 'birth';
        else if (rand < 0.67) return 'death';
        else return 'jump';
    }
}

// Run RJMCMC step
function runRJMCMCStep() {
    // Run RJMCMC step for each chain
    for (let chainId = 0; chainId < rjmcState.chains.length; chainId++) {
        const chain = rjmcState.chains[chainId];
        const currentK = chain.currentJump.p.length;
        const proposalType = sampleProposalType(currentK);
        
        let proposedJump = null;
        let logAcceptanceRatio = 0;
        
        if (proposalType === 'birth' && currentK < 8) {
            proposedJump = sampleBirthProposal(chain.currentJump);
            logAcceptanceRatio = Math.log(1 / (currentK * (1/20))); // Birth proposal probability
        } else if (proposalType === 'death' && currentK > 2) {
            const indexToRemove = Math.floor(Math.random() * currentK);
            proposedJump = sampleDeathProposal(chain.currentJump, indexToRemove);
            logAcceptanceRatio = Math.log(currentK * (1/20)); // Death proposal probability
        } else {
            // Jump move
            const indexToUpdate = Math.floor(Math.random() * currentK);
            proposedJump = sampleJumpProposal(chain.currentJump, indexToUpdate);
            logAcceptanceRatio = 0; // Symmetric proposal
        }
        
        if (proposedJump) {
            const proposedLogPost = evaluateLogPosterior(chain.currentParams, proposedJump);
            const acceptanceRatio = Math.exp(proposedLogPost - chain.logPosterior + logAcceptanceRatio);
            
            if (Math.random() < acceptanceRatio) {
                // Accept proposal
                chain.currentJump = proposedJump;
                chain.logPosterior = proposedLogPost;
                chain.acceptedMoves++;
            }
        }
        
        // Update chain history
        chain.currentK = chain.currentJump.p.length;
        chain.kHistory.push(chain.currentK);
        chain.logPostHistory.push(chain.logPosterior);
        
        // Update kProbabilities immediately for accurate tab ordering
        // This is just for tracking, the real calculation happens in updateComponentDistribution
        rjmcState.kProbabilities[chain.currentK] = (rjmcState.kProbabilities[chain.currentK] || 0) + 1;
        
        // Limit history length to prevent memory issues
        if (chain.kHistory.length > rjmcState.maxHistoryLength) {
            chain.kHistory.shift();
            chain.logPostHistory.shift();
        }
        
        // Store parameter traces for all components (optimized)
        if (chain.currentJump.mu && chain.currentJump.mu.length > 0) {
            // Store references instead of deep copies for performance
            chain.muHistory.push(chain.currentJump.mu);
            chain.sigmaHistory.push(chain.currentJump.sigma);
            
            // Limit history length to prevent memory issues
            if (chain.muHistory.length > rjmcState.maxHistoryLength) {
                chain.muHistory.shift(); // Remove oldest entry
                chain.sigmaHistory.shift();
            }
        }
    }
    
    // Update combined state by aggregating all chains
    rjmcState.iterations++;
    
    // Periodic cleanup to prevent memory buildup
    if (rjmcState.iterations % rjmcState.cleanupFrequency === 0) {
        cleanupOldData();
    }
    
    // kProbabilities are now updated every iteration in runRJMCMCStep
    
    // Use first chain for current state display (representative)
    const firstChain = rjmcState.chains[0];
    rjmcState.currentJump = firstChain.currentJump;
    rjmcState.logPosterior = firstChain.logPosterior;
    rjmcState.currentK = firstChain.currentK;
    rjmcState.kHistory.push(firstChain.currentK);
    rjmcState.logPostHistory.push(firstChain.logPosterior);
    
    // Store best fit for current K if it's better than previous
    const currentK = rjmcState.currentK;
    if (!rjmcState.topKFits[currentK] || rjmcState.logPosterior > rjmcState.topKFits[currentK].logPosterior) {
        rjmcState.topKFits[currentK] = {
            jump: JSON.parse(JSON.stringify(rjmcState.currentJump)), // Deep copy
            logPosterior: rjmcState.logPosterior
        };
    }
    
    // Store all posterior samples for each K from all chains (after burn-in)
    const burninInput = document.getElementById('burnin');
    const burnin = burninInput ? parseInt(burninInput.value) || 100 : 100;
    
    if (rjmcState.iterations > burnin) {
        rjmcState.chains.forEach(chain => {
            const k = chain.currentK;
            if (!rjmcState.kPosteriorSamples[k]) {
                rjmcState.kPosteriorSamples[k] = [];
            }
            rjmcState.kPosteriorSamples[k].push({
                jump: JSON.parse(JSON.stringify(chain.currentJump)),
                logPosterior: chain.logPosterior
            });
        });
        
        // Debug logging
        if (rjmcState.iterations % 1000 === 0) {
            console.log('[RJMCMC] kPosteriorSamples status:', Object.keys(rjmcState.kPosteriorSamples).map(k => 
                `${k}: ${rjmcState.kPosteriorSamples[k].length} samples`
            ));
        }
    }
}

// Run RJMCMC analysis with streaming updates
async function runRJMCMCAnalysis(settings = {}) {
    if (isAnalysisRunning) {
        console.log('[RJMCMC] Analysis already running, skipping...');
        return;
    }
    
    if (!currentData || !currentData.obs || currentData.obs.length === 0) {
        throw new Error('No data available for analysis');
    }
    
    isAnalysisRunning = true;
    
    try {
        console.log('[RJMCMC] Starting RJMCMC analysis...');
        
        // Default settings
        const defaultSettings = {
            iterations: 20000,
            burnin: 10000,
            updateInterval: 50
        };
        
        const finalSettings = { ...defaultSettings, ...settings };
        
        // Initialize RJMCMC with multiple chains
        initializeRJMCMC();
        
        // Show initial progress
        updateProgress(0, finalSettings.iterations, 'Initializing RJMCMC...');
        showProgressBar();
        updateProgressBar(0, finalSettings.iterations);
        
        // Run RJMCMC iterations
        console.log(`[RJMCMC] Starting ${finalSettings.iterations} iterations`);
        for (let i = 0; i < finalSettings.iterations && isAnalysisRunning; i++) {
            try {
            runRJMCMCStep();
            } catch (stepError) {
                console.error(`[RJMCMC] Error in step ${i}:`, stepError);
                throw stepError;
            }
            
            // Debug logging for iteration count
            if (i % 1000 === 0 && i > 0) {
                console.log(`[RJMCMC] Completed ${i} iterations, isAnalysisRunning: ${isAnalysisRunning}`);
            }
            
            // Update progress and visualizations (optimized frequency)
            if (i % finalSettings.updateInterval === 0 || i < 100) {
                updateProgress(i + 1, finalSettings.iterations, `Running RJMCMC... K=${rjmcState.currentK}`);
                updateProgressBar(i + 1, finalSettings.iterations);
                
                // Update visualizations with fixed frequency for better performance
                const updateFreq = i < 100 ? 10 : 100; // Every 10 iterations for first 100, then every 100
                if (i % updateFreq === 0 || i < 100 || i === finalSettings.iterations - 1) {
                    updateVisualizations(); // Update all visualizations
                }
                
                // Always add a small delay for UI responsiveness and to prevent browser freezing
                await new Promise(resolve => setTimeout(resolve, 1));
            }
            
            // Update convergence diagnostics every 100 iterations (independent of visualization updates)
            if (i % 100 === 0 && i > 0) {
                console.log(`[RJMCMC] Updating convergence diagnostics at iteration ${i}`);
                updateConvergenceDiagnostics();
            }
            
            // Update mixture fit and component distribution every 500 iterations for better performance
            if (i % 500 === 0 && i > 0) {
                updateMixtureFit();
                updateComponentDistribution();
                
                // Show posterior distribution container when we have data
                const posteriorContainer = document.getElementById('posteriorDistContainer');
                if (posteriorContainer && rjmcState.kHistory.length > 0) {
                    posteriorContainer.style.display = 'block';
                }
            }
        }
        
        if (isAnalysisRunning) {
            console.log(`[RJMCMC] Analysis completed successfully after ${finalSettings.iterations} iterations`);
            updateProgress(finalSettings.iterations, finalSettings.iterations, 'Analysis complete!');
            updateProgressBar(finalSettings.iterations, finalSettings.iterations);
            
            // Update convergence diagnostics one final time
            console.log('[RJMCMC] Final convergence diagnostics update');
            updateConvergenceDiagnostics();
            
        // Set MCMC as not running and update parameter trace plots
        setMCMCRunningState(false);
        updateParameterTracePlots(); // Update parameter plots once at the end
        
        // Update mixture fit and component distribution
        updateMixtureFit();
        updateComponentDistribution();
        
        // Update mixture tabs to show available K values
        updateMixtureTabs();
        
        // Force show posterior distribution container
        const posteriorContainer = document.getElementById('posteriorDistContainer');
        if (posteriorContainer) {
            posteriorContainer.style.display = 'block';
            console.log('[RJMCMC] Posterior distribution container made visible');
        }
            
        updateFinalResults();
        
        // Hide progress bar after completion
        setTimeout(() => {
            hideProgressBar();
        }, 2000); // Hide after 2 seconds to show completion
        
        } else {
            console.log(`[RJMCMC] Analysis was stopped early. isAnalysisRunning: ${isAnalysisRunning}`);
        }
        
    } catch (error) {
        console.error('[RJMCMC] Analysis failed:', error);
        console.error('[RJMCMC] Error stack:', error.stack);
        updateProgress(0, 0, 'Analysis failed: ' + error.message);
        hideProgressBar(); // Hide progress bar on error
        throw error;
        
    } finally {
        isAnalysisRunning = false;
        cleanupAutoTabSwitching();
    }
}

// Update progress display
function updateProgress(current, total, message) {
    const progressElement = document.getElementById('mcmc-progress');
    const progressBar = document.getElementById('progress-bar');
    
    const percentage = total > 0 ? (current / total * 100).toFixed(1) : 0;
    
    if (progressElement && progressBar) {
        progressElement.textContent = `Step ${current} / ${total} (${percentage}%) - ${message}`;
        progressBar.style.width = percentage + '%';
    }
    
    console.log(`[PROGRESS] ${current}/${total} (${percentage}%) - ${message}`);
}

// Update visualizations
function updateVisualizations() {
    try {
        updateComponentTracePlot();
        updateLogPosteriorPlot();
        updateComponentDistribution();
        updateMixtureFit();
    } catch (error) {
        console.error('[RJMCMC] Error in updateVisualizations:', error);
    }
}

// Update final results
function updateFinalResults() {
    try {
        console.log('[RJMCMC] updateFinalResults called');
        console.log('[RJMCMC] kProbabilities:', rjmcState.kProbabilities);
        console.log('[RJMCMC] kHistory length:', rjmcState.kHistory.length);
        
        // Calculate final K probabilities
        const totalIterations = rjmcState.kHistory.length;
        const kProbs = {};
        for (const [k, count] of Object.entries(rjmcState.kProbabilities)) {
            kProbs[k] = count / totalIterations;
        }
        
        console.log('[RJMCMC] Calculated kProbs:', kProbs);
        
        // Find most probable K
        const mostProbableK = Object.keys(kProbs).reduce((a, b) => 
            kProbs[a] > kProbs[b] ? a : b
        );
        
        console.log('[RJMCMC] Most probable K:', mostProbableK);
        
        // Update displays
        const kElement = document.getElementById('current-k');
        if (kElement) {
            kElement.textContent = mostProbableK;
        }
        
        updateKProbabilities(kProbs);
        
        // Force update all charts
        console.log('[RJMCMC] Updating all charts...');
        updateComponentTracePlot();
        updateLogPosteriorPlot();
        updateComponentDistribution();
        updateMixtureFit();
        updateParameterTracePlots();
        
        // Force update mixture tabs
        updateMixtureTabs();
        
        console.log('[RJMCMC] Final results updated successfully');
        
    } catch (error) {
        console.error('[RJMCMC] Failed to update final results:', error);
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

// Update all charts
function updateAllCharts() {
    updateComponentTracePlot();
    updateLogPosteriorPlot();
    updateComponentDistribution();
    updateMixtureFit();
}

// Test function to create a simple plot
function testComponentTracePlot() {
    console.log('[RJMCMC] Testing component trace plot...');
    
    const testData = [{
        x: [1, 2, 3, 4, 5],
        y: [2, 3, 2, 4, 3],
        type: 'scatter',
        mode: 'lines',
        name: 'Test Chain',
        line: { color: '#000000', width: 2 }
    }];
    
    const layout = {
        title: 'Test Plot',
        xaxis: { title: 'Iteration', autorange: true },
        yaxis: { title: 'K', autorange: true, fixedrange: false }
    };
    
    try {
        Plotly.newPlot('componentTraceChart', testData, layout);
        console.log('[RJMCMC] Test plot created successfully');
    } catch (error) {
        console.error('[RJMCMC] Test plot failed:', error);
    }
}

// Update component trace plot
function updateComponentTracePlot() {
    console.log('[RJMCMC] updateComponentTracePlot called');
    console.log('[RJMCMC] rjmcState.kHistory:', rjmcState.kHistory);
    console.log('[RJMCMC] rjmcState.chains:', rjmcState.chains);
    
    if (!rjmcState.kHistory || rjmcState.kHistory.length === 0) {
        console.log('[RJMCMC] No kHistory data, showing test plot instead');
        testComponentTracePlot();
        return;
    }
    
    try {
        // Create Plotly traces for all chains
        const plotData = [];
        const colors = ['#000000', '#666666', '#999999', '#CCCCCC'];
        
        for (let chainId = 0; chainId < rjmcState.chains.length; chainId++) {
            const chain = rjmcState.chains[chainId];
            if (chain.kHistory.length === 0) continue;
            
            // Use smart data range based on MCMC state
            const { data: smoothedK, labels: smoothedLabels } = getDataRange(
                chain.kHistory, 
                rjmcState.isRunning, 
                rjmcState.showBurnin
            );
            
            plotData.push({
                x: smoothedLabels,
                y: smoothedK,
                type: 'scatter',
                mode: 'lines',
                name: `Chain ${chainId + 1}`,
                line: {
                    color: colors[chainId],
                    width: 2
                },
                marker: {
                    size: 3
                }
            });
        }
        
        console.log('[RJMCMC] plotData created:', plotData);
        console.log('[RJMCMC] Number of traces:', plotData.length);
        
        const layout = {
            title: {
                text: 'Number of Components',
                font: { size: 14, family: 'Avenir, sans-serif' }
            },
            xaxis: {
                title: 'Iteration',
                showgrid: true,
                gridcolor: '#e0e0e0',
                font: { family: 'Avenir, sans-serif' },
                autorange: true,
                fixedrange: false
            },
            yaxis: {
                title: 'Number of Components',
                showgrid: true,
                gridcolor: '#e0e0e0',
                font: { family: 'Avenir, sans-serif' },
                autorange: true,
                fixedrange: false,
                dtick: 1
            },
            margin: { l: 60, r: 20, t: 40, b: 50 },
            plot_bgcolor: 'rgba(0,0,0,0)',
            paper_bgcolor: 'rgba(0,0,0,0)',
            font: { family: 'Avenir, sans-serif' },
            legend: {
                font: { family: 'Avenir, sans-serif' }
            }
        };
        
        const config = {
            responsive: true,
            displayModeBar: false
        };
        
        const element = document.getElementById('componentTraceChart');
        console.log('[RJMCMC] Canvas element found:', element);
        
        if (!element) {
            console.error('[RJMCMC] Canvas element componentTraceChart not found!');
            return;
        }
        
        try {
            Plotly.newPlot('componentTraceChart', plotData, layout, config);
            console.log('[RJMCMC] Component trace plot created successfully');
        } catch (error) {
            console.error('[RJMCMC] Failed to create component trace Plotly chart:', error);
        }
        
    } catch (error) {
        console.error('[RJMCMC] Failed to update component trace plot:', error);
    }
}

// Update log posterior plot
function updateLogPosteriorPlot() {
    if (!rjmcState.logPostHistory || rjmcState.logPostHistory.length === 0) return;
    
    try {
        // Create Plotly traces for all chains
        const plotData = [];
        const colors = ['#000000', '#666666', '#999999', '#CCCCCC'];
        
        for (let chainId = 0; chainId < rjmcState.chains.length; chainId++) {
            const chain = rjmcState.chains[chainId];
            if (chain.logPostHistory.length === 0) continue;
            
            // Use smart data range based on MCMC state
            const { data: smoothedLogPost, labels: smoothedLabels } = getDataRange(
                chain.logPostHistory, 
                rjmcState.isRunning, 
                rjmcState.showBurnin
            );
            
            plotData.push({
                x: smoothedLabels,
                y: smoothedLogPost,
                type: 'scatter',
                mode: 'lines',
                name: `Chain ${chainId + 1}`,
                line: {
                    color: colors[chainId],
                    width: 2
                },
                marker: {
                    size: 3
                }
            });
        }
        
        const layout = {
            title: {
                text: 'Log Posterior',
                font: { size: 14, family: 'Avenir, sans-serif' }
            },
            xaxis: {
                title: 'Iteration',
                showgrid: true,
                gridcolor: '#e0e0e0',
                font: { family: 'Avenir, sans-serif' },
                autorange: true,
                fixedrange: false
            },
            yaxis: {
                title: 'Log Posterior',
                showgrid: true,
                gridcolor: '#e0e0e0',
                font: { family: 'Avenir, sans-serif' },
                autorange: true,
                fixedrange: false
            },
            margin: { l: 60, r: 20, t: 40, b: 50 },
            plot_bgcolor: 'rgba(0,0,0,0)',
            paper_bgcolor: 'rgba(0,0,0,0)',
            font: { family: 'Avenir, sans-serif' },
            legend: {
                font: { family: 'Avenir, sans-serif' }
            }
        };
        
        const config = {
            responsive: true,
            displayModeBar: false
        };
        
        const element = document.getElementById('logPosteriorChart');
        if (!element) {
            console.error('[RJMCMC] Canvas element logPosteriorChart not found!');
            return;
        }
        
        try {
            Plotly.newPlot('logPosteriorChart', plotData, layout, config);
            console.log('[RJMCMC] Log posterior plot created successfully');
        } catch (error) {
            console.error('[RJMCMC] Failed to create log posterior Plotly chart:', error);
        }
        
    } catch (error) {
        console.error('[RJMCMC] Failed to update log posterior plot:', error);
    }
}

// Update component distribution
function updateComponentDistribution() {
    console.log('[RJMCMC] updateComponentDistribution called');
    
    const chartElement = document.getElementById('componentDistChart');
    if (!chartElement) {
        console.error('[RJMCMC] componentDistChart element not found');
        return;
    }
    
    const ctx = chartElement.getContext('2d');
    
    if (componentDistChart) {
        componentDistChart.destroy();
        componentDistChart = null;
    }
    
    // Collect all K values from all chains
    let allKValues = [];
    
    if (rjmcState.chains && rjmcState.chains.length > 0) {
        console.log('[RJMCMC] Collecting K values from all chains');
        rjmcState.chains.forEach((chain, chainIndex) => {
            if (chain.kHistory && chain.kHistory.length > 0) {
                allKValues = allKValues.concat(chain.kHistory);
                console.log(`[RJMCMC] Chain ${chainIndex + 1}: ${chain.kHistory.length} K values`);
            }
        });
    } else if (rjmcState.kHistory && rjmcState.kHistory.length > 0) {
        console.log('[RJMCMC] Using single chain kHistory');
        allKValues = rjmcState.kHistory;
    }
    
    console.log('[RJMCMC] Total K values collected:', allKValues.length);
    console.log('[RJMCMC] K values sample:', allKValues.slice(0, 10));
    
    if (allKValues.length === 0) {
        console.log('[RJMCMC] No K values available, creating empty chart');
        componentDistChart = new Chart(ctx, {
            type: 'bar',
            data: {
                labels: ['No data yet'],
                datasets: [{
                    label: 'Frequency',
                    data: [0],
                    backgroundColor: 'rgba(0,0,0,0.1)',
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
                        max: 100,
                        title: {
                            display: true,
                            text: 'Frequency (%)'
                        }
                    },
                    x: {
                        title: {
                            display: true,
                            text: 'Number of Components (K)'
                        }
                    }
                },
                plugins: {
                    title: {
                        display: true,
                        text: 'Posterior Distribution of Components'
                    }
                }
            }
        });
        return;
    }
    
    // Count frequencies of each K value
    const kCounts = {};
    allKValues.forEach(k => {
        kCounts[k] = (kCounts[k] || 0) + 1;
    });
    
    console.log('[RJMCMC] K counts:', kCounts);
    
    // Get sorted K values
    const kValues = Object.keys(kCounts).map(Number).sort((a, b) => a - b);
    const totalCount = allKValues.length;
    const frequencies = kValues.map(k => (kCounts[k] / totalCount) * 100);
    
    console.log('[RJMCMC] Creating histogram with:', {
        kValues: kValues,
        frequencies: frequencies,
        totalCount: totalCount
    });
    
    // Check if Chart.js is available
    if (typeof Chart === 'undefined') {
        console.error('[RJMCMC] Chart.js is not loaded!');
        return;
    }
    
    try {
        componentDistChart = new Chart(ctx, {
            type: 'bar',
            data: {
                labels: kValues.map(k => `K=${k}`),
                datasets: [{
                    label: 'Frequency (%)',
                    data: frequencies,
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
                        max: Math.max(...frequencies) * 1.1,
                        title: {
                            display: true,
                            text: 'Frequency (%)'
                        }
                    },
                    x: {
                        title: {
                            display: true,
                            text: 'Number of Components (K)'
                        }
                    }
                },
                plugins: {
                    title: {
                        display: true,
                        text: 'Posterior Distribution of Components'
                    }
                }
            }
        });
        
        console.log('[RJMCMC] Chart created successfully');
        
        // Show the posterior distribution container
        const container = document.getElementById('posteriorDistContainer');
        if (container) {
            container.style.display = 'block';
        }
        
    } catch (error) {
        console.error('[RJMCMC] Error creating chart:', error);
    }
}

// Update tab labels with top 3 K values
function updateTabLabels(kValues, probabilities) {
    // Find top 3 K values by probability
    const kWithProbs = kValues.map((k, i) => ({ k, prob: probabilities[i] }))
        .sort((a, b) => b.prob - a.prob)
        .slice(0, 3);
    
    // Update tab buttons
    const tabButtons = ['top1', 'top2', 'top3'];
    tabButtons.forEach((tabId, index) => {
        const button = document.getElementById(`tab-${tabId}`);
        if (button && kWithProbs[index]) {
            const k = kWithProbs[index].k;
            const prob = kWithProbs[index].prob.toFixed(1);
            button.textContent = `K=${k} (${prob}%)`;
        }
    });
}

// Calculate posterior means for each K
function calculatePosteriorMeans() {
    for (const k in rjmcState.kPosteriorSamples) {
        const samples = rjmcState.kPosteriorSamples[k];
        if (samples.length === 0) continue;
        
        const numComponents = parseInt(k);
        const posteriorMeans = {
            p: new Array(numComponents).fill(0),
            mu: new Array(numComponents).fill(0),
            sigma: new Array(numComponents).fill(0)
        };
        
        // Calculate means across all samples, ordering each sample's means first
        samples.forEach(sample => {
            // Order the means from lowest to highest for this sample
            const muValues = sample.jump.mu.slice();
            const sigmaValues = sample.jump.sigma.slice();
            const pValues = sample.jump.p.slice();
            
            // Create indices sorted by mean values
            const sortedIndices = muValues
                .map((val, idx) => ({ val, idx }))
                .sort((a, b) => a.val - b.val)
                .map(item => item.idx);
            
            // Add ordered values to posterior means
            for (let i = 0; i < numComponents; i++) {
                const originalIdx = sortedIndices[i];
                posteriorMeans.p[i] += pValues[originalIdx];
                posteriorMeans.mu[i] += muValues[originalIdx];
                posteriorMeans.sigma[i] += sigmaValues[originalIdx];
            }
        });
        
        // Normalize by number of samples
        const numSamples = samples.length;
        for (let i = 0; i < numComponents; i++) {
            posteriorMeans.p[i] /= numSamples;
            posteriorMeans.mu[i] /= numSamples;
            posteriorMeans.sigma[i] /= numSamples;
        }
        
        rjmcState.kPosteriorMeans[k] = posteriorMeans;
    }
}

// Calculate Effective Sample Size
function calculateESS(chainData) {
    if (chainData.length < 10) return 0;
    
    const n = chainData.length;
    const mean = chainData.reduce((sum, val) => sum + val, 0) / n;
    
    // Calculate autocorrelation
    let autocorr = [];
    for (let lag = 0; lag < Math.min(n/4, 100); lag++) {
        let numerator = 0;
        let denominator = 0;
        
        for (let i = 0; i < n - lag; i++) {
            numerator += (chainData[i] - mean) * (chainData[i + lag] - mean);
        }
        
        for (let i = 0; i < n; i++) {
            denominator += (chainData[i] - mean) * (chainData[i] - mean);
        }
        
        if (denominator === 0) break;
        const rho = numerator / denominator;
        autocorr.push(rho);
        
        // Stop if autocorrelation becomes negative
        if (rho < 0) break;
    }
    
    // Calculate ESS
    let sumRho = 0;
    for (let i = 1; i < autocorr.length; i++) {
        sumRho += autocorr[i];
    }
    
    const ess = n / (1 + 2 * sumRho);
    return Math.max(0, Math.round(ess));
}

// Calculate Gelman-Rubin statistic
function calculateGelmanRubin(chains) {
    if (chains.length < 2) return 1.0;
    
    const M = chains.length;
    const n = chains[0].length;
    
    if (n < 10) return 1.0;
    
    // Calculate within-chain means and variances
    const chainMeans = chains.map(chain => {
        return chain.reduce((sum, val) => sum + val, 0) / n;
    });
    
    const chainVars = chains.map(chain => {
        const mean = chainMeans[chains.indexOf(chain)];
        return chain.reduce((sum, val) => sum + Math.pow(val - mean, 2), 0) / (n - 1);
    });
    
    // Calculate overall mean
    const overallMean = chainMeans.reduce((sum, mean) => sum + mean, 0) / M;
    
    // Calculate between-chain variance
    const B = (n / (M - 1)) * chainMeans.reduce((sum, mean) => sum + Math.pow(mean - overallMean, 2), 0);
    
    // Calculate within-chain variance
    const W = chainVars.reduce((sum, var_) => sum + var_, 0) / M;
    
    // Calculate pooled variance
    const V = ((n - 1) / n) * W + (1 / n) * B;
    
    // Calculate R-hat
    const Rhat = Math.sqrt(V / W);
    
    return Rhat;
}

// Update convergence diagnostics
function updateConvergenceDiagnostics() {
    console.log('[RJMCMC] updateConvergenceDiagnostics called');
    if (rjmcState.chains.length === 0) {
        console.log('[RJMCMC] No chains available for convergence diagnostics');
        return;
    }
    
    // Get most common K
    console.log('[RJMCMC] kProbabilities:', rjmcState.kProbabilities);
    const mostCommonK = getMostCommonK();
    console.log('[RJMCMC] Most common K:', mostCommonK);
    if (!mostCommonK) {
        console.log('[RJMCMC] No most common K found, skipping convergence diagnostics');
        return;
    }
    
    // Generate parameter traces using the shared function
    const traces = generateParameterTraces();
    if (!traces) return;
    
    const { meanTraces, sigmaTraces } = traces;
    
    // Calculate ESS for K using last 5000 steps (excluding burn-in)
    const combinedKHistory = [];
    const maxSteps = 5000; // Use last 5000 steps for convergence diagnostics
    const burnin = 100; // Burn-in period
    
    rjmcState.chains.forEach(chain => {
        // Remove burn-in and take last 5000 steps
        const postBurnin = chain.kHistory.slice(burnin);
        const lastSteps = postBurnin.slice(-maxSteps);
        combinedKHistory.push(...lastSteps);
    });
    const kESS = calculateESS(combinedKHistory);
    
    // Calculate Gelman-Rubin for K using last 5000 steps (excluding burn-in)
    const kTraces = [];
    rjmcState.chains.forEach(chain => {
        const postBurnin = chain.kHistory.slice(burnin);
        const lastSteps = postBurnin.slice(-maxSteps);
        if (lastSteps.length > 0) {
            kTraces.push(lastSteps);
        }
    });
    const kRhat = kTraces.length > 1 ? calculateGelmanRubin(kTraces) : 1.0;
    
    // Debug logging for convergence diagnostics
    console.log(`[RJMCMC] Convergence diagnostics: Using ${combinedKHistory.length} total steps from ${kTraces.length} chains (last ${maxSteps} steps excluding burn-in)`);
    
    // Calculate acceptance rates for last 5000 steps (excluding burn-in)
    let totalMoves = 0;
    let totalAccepted = 0;
    
    rjmcState.chains.forEach(chain => {
        // Calculate acceptance rate for last 5000 steps
        const postBurnin = chain.kHistory.slice(burnin);
        const lastSteps = postBurnin.slice(-maxSteps);
        const movesInWindow = lastSteps.length;
        const acceptedInWindow = Math.min(chain.acceptedMoves, movesInWindow); // Approximate accepted moves in window
        
        totalMoves += movesInWindow;
        totalAccepted += acceptedInWindow;
    });
    
    const acceptanceRate = totalMoves > 0 ? (totalAccepted / totalMoves * 100).toFixed(1) : 0;
    
    // Update displays
    console.log(`[RJMCMC] Updating displays: ESS=${kESS}, Rhat=${kRhat}, Acceptance=${acceptanceRate}%`);
    updateESSDisplay(kESS);
    updateGelmanRubinDisplay(kRhat);
    updateAcceptanceRatesDisplay(acceptanceRate);
    
    // Update parameter trace plots (commented out for performance - too expensive for real-time)
    // updateParameterTracePlots(meanTraces, sigmaTraces, mostCommonK);
}

// Update ESS display
function updateESSDisplay(ess) {
    console.log(`[RJMCMC] updateESSDisplay called with ESS: ${ess}`);
    
    // Update the main ESS display
    const essValue = document.getElementById('essValue');
    if (essValue) {
        essValue.textContent = ess;
    }
    
    // Also update the detailed display if it exists
    const essDisplay = document.getElementById('ess-display');
    if (essDisplay) {
        let status = 'convergence-poor';
        if (ess > 1000) status = 'convergence-good';
        else if (ess > 100) status = 'convergence-warning';
        
        essDisplay.innerHTML = `
            <div class="diagnostic-value ${status}">${ess}</div>
            <div class="diagnostic-label">K (Number of Components)</div>
            <small>ESS > 1000: Good, ESS > 100: Acceptable, ESS < 100: Poor</small>
        `;
    }
}

// Update Gelman-Rubin display
function updateGelmanRubinDisplay(rhat) {
    // Update the main R-hat display
    const rhatValue = document.getElementById('rhatValue');
    const rhatLabel = document.getElementById('rhatLabel');
    if (rhatValue) {
        rhatValue.textContent = rhat.toFixed(3);
    }
    if (rhatLabel) {
        if (rhat < 1.1) {
            rhatLabel.textContent = '✓ Converged (R̂ < 1.1)';
            rhatLabel.style.color = '#28a745';
        } else if (rhat < 1.2) {
            rhatLabel.textContent = '⚠ Acceptable (R̂ < 1.2)';
            rhatLabel.style.color = '#ffc107';
        } else {
            rhatLabel.textContent = '❌ Poor convergence (R̂ > 1.2)';
            rhatLabel.style.color = '#dc3545';
        }
    }
    
    // Also update the detailed display if it exists
    const grDisplay = document.getElementById('gelman-rubin-display');
    if (grDisplay) {
        let status = 'convergence-poor';
        let convergenceMessage = '';
        if (rhat < 1.1) {
            status = 'convergence-good';
            convergenceMessage = '<div style="color: #28a745; font-weight: bold; margin-top: 5px;">✓ Evidence of convergence!</div>';
        } else if (rhat < 1.2) {
            status = 'convergence-warning';
        }
        
        grDisplay.innerHTML = `
            <div class="diagnostic-value ${status}">${rhat.toFixed(3)}</div>
            <div class="diagnostic-label">K (Number of Components)</div>
            <small>R̂ < 1.1: Good, R̂ < 1.2: Acceptable, R̂ > 1.2: Poor</small>
            ${convergenceMessage}
        `;
    }
}

// Update acceptance rates display
function updateAcceptanceRatesDisplay(rate) {
    // Update the main acceptance rate display
    const acceptanceValue = document.getElementById('acceptanceValue');
    if (acceptanceValue) {
        acceptanceValue.textContent = rate + '%';
    }
    
    // Also update the detailed display if it exists
    const arDisplay = document.getElementById('acceptance-rates-display');
    if (arDisplay) {
        let status = 'convergence-poor';
        if (rate > 20 && rate < 50) status = 'convergence-good';
        else if (rate > 10 && rate < 70) status = 'convergence-warning';
        
        arDisplay.innerHTML = `
            <div class="diagnostic-value ${status}">${rate}%</div>
            <div class="diagnostic-label">Overall Acceptance Rate</div>
            <small>20-50%: Good, 10-70%: Acceptable, Outside: Poor</small>
        `;
    }
}

// Generate parameter traces from current state
function generateParameterTraces() {
    if (rjmcState.chains.length === 0) return null;
    
    // Get most common K
    const mostCommonK = getMostCommonK();
    if (!mostCommonK) return null;
    
    // Collect parameter traces from all chains for the most common K
    const meanTraces = [];
    const sigmaTraces = [];
    const kTraces = [];
    
    rjmcState.chains.forEach(chain => {
        // Remove burn-in (first 100 iterations)
        const burnin = Math.min(100, chain.kHistory.length);
        const kHistory = chain.kHistory.slice(burnin);
        const meanHistory = chain.muHistory.slice(burnin);
        const sigmaHistory = chain.sigmaHistory.slice(burnin);
        
        if (kHistory.length > 0) {
            kTraces.push(kHistory);
        }
        
        // Extract traces for the most common K only
        const filteredMeanHistory = [];
        const filteredSigmaHistory = [];
        
        for (let i = 0; i < meanHistory.length; i++) {
            if (i < kHistory.length && kHistory[i] === mostCommonK) {
                // Get the component values for this K
                const muValues = meanHistory[i];
                const sigmaValues = sigmaHistory[i];
                
                if (muValues && sigmaValues && muValues.length === mostCommonK) {
                    // Rank the means in ascending order (mean1 < mean2 < mean3, etc.)
                    const rankedIndices = muValues
                        .map((val, idx) => ({ val, idx }))
                        .sort((a, b) => a.val - b.val)
                        .map(item => item.idx);
                    
                    // Reorder both means and sigmas according to the ranking
                    const rankedMuValues = rankedIndices.map(idx => muValues[idx]);
                    const rankedSigmaValues = rankedIndices.map(idx => sigmaValues[idx]);
                    
                    filteredMeanHistory.push(rankedMuValues);
                    filteredSigmaHistory.push(rankedSigmaValues);
                }
            }
        }
        
        if (filteredMeanHistory.length > 0) {
            meanTraces.push(filteredMeanHistory);
            sigmaTraces.push(filteredSigmaHistory);
        }
    });
    
    return {
        meanTraces,
        sigmaTraces,
        mostCommonK
    };
}

// Update parameter trace plots (now faceted)
function updateParameterTracePlots(meanTraces, sigmaTraces, mostCommonK) {
    // If no parameters provided, generate them from current state
    if (!meanTraces || !sigmaTraces || !mostCommonK) {
        const traces = generateParameterTraces();
        if (traces) {
            meanTraces = traces.meanTraces;
            sigmaTraces = traces.sigmaTraces;
            mostCommonK = traces.mostCommonK;
        } else {
            console.log('[RJMCMC] No parameter trace data available');
            return;
        }
    }
    updateFacetedParameterPlots(meanTraces, sigmaTraces, mostCommonK);
}

// Update faceted parameter plots
function updateFacetedParameterPlots(meanTraces, sigmaTraces, mostCommonK) {
    if (meanTraces.length === 0 || !mostCommonK) return;
    
    try {
        const container = document.getElementById('parameterTraceGrid');
        if (!container) {
            console.log('[RJMCMC] parameterTraceGrid container not found');
            return;
        }
        
        // Clear existing plots
        container.innerHTML = '';
        
        // Set up container layout - use vertical stacking (rows)
        container.style.display = 'flex';
        container.style.flexDirection = 'column';
        container.style.gap = '20px';
        container.style.padding = '10px 0';
        
        // Create faceted plots for each component (ranked by mean value)
        for (let compId = 0; compId < mostCommonK; compId++) {
            // Create container for this component
            const componentContainer = document.createElement('div');
            componentContainer.className = 'parameter-trace-item';
            
            componentContainer.innerHTML = `
                <h4>Component ${compId + 1}</h4>
                <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 15px; margin-top: 10px; width: 100%;">
                    <div style="width: 100%;">
                        <h5 style="font-size: 0.9em; margin-bottom: 8px; text-align: center; font-weight: 600;">Mean (μ${compId + 1})</h5>
                        <canvas id="meanChart_${compId}" class="parameter-trace-chart" style="width: 100%; height: 200px;"></canvas>
                    </div>
                    <div style="width: 100%;">
                        <h5 style="font-size: 0.9em; margin-bottom: 8px; text-align: center; font-weight: 600;">Std Dev (σ${compId + 1})</h5>
                        <canvas id="sigmaChart_${compId}" class="parameter-trace-chart" style="width: 100%; height: 200px;"></canvas>
                    </div>
                </div>
            `;
            container.appendChild(componentContainer);
            
            // Create mean chart for this component
            createComponentChart(`meanChart_${compId}`, meanTraces, compId, 'Component Mean', '#FF6B6B');
            
            // Create sigma chart for this component
            createComponentChart(`sigmaChart_${compId}`, sigmaTraces, compId, 'Component Sigma', '#4ECDC4');
        }
        
    } catch (error) {
        console.error('[RJMCMC] Failed to update faceted parameter plots:', error);
    }
}

// Create individual component chart
function createComponentChart(canvasId, traces, compId, title, color) {
    try {
        const ctx = document.getElementById(canvasId).getContext('2d');
        
        const chainColors = ['#000000', '#666666', '#999999', '#CCCCCC'];
        const datasets = [];
        
        traces.forEach((chainTrace, chainId) => {
            if (chainTrace.length === 0) return;
            
            // Extract the specific component's trace
            const componentTrace = chainTrace.map(iteration => {
                return iteration && iteration[compId] !== undefined ? iteration[compId] : null;
            }).filter(val => val !== null);
            
            // Use smart data range based on MCMC state
            const { data: displayTrace, labels: traceLabels } = getDataRange(
                componentTrace, 
                rjmcState.isRunning, 
                rjmcState.showBurnin
            );
            
            if (displayTrace.length > 0) {
                datasets.push({
                    label: `Chain ${chainId + 1}`,
                    data: displayTrace,
                    borderColor: chainColors[chainId % chainColors.length], // Use modulo to handle more than 4 chains
                    backgroundColor: chainColors[chainId % chainColors.length] + '20',
                    borderWidth: 1,
                    fill: false,
                    pointRadius: 0,
                    pointHoverRadius: 2,
                    borderDash: chainId > 0 ? [3, 3] : []
                });
            }
        });
        
        // Create proper iteration labels for the x-axis
        const iterationLabels = datasets.length > 0 ? 
            datasets[0].data.map((_, i) => i) : []; // Use the labels from getDataRange
        
        new Chart(ctx, {
            type: 'line',
            data: {
                labels: iterationLabels,
                datasets: datasets
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                plugins: {
                    legend: {
                        display: false // Hide legend for cleaner look
                    }
                },
                scales: {
                    y: {
                        beginAtZero: false,
                        title: {
                            display: false
                        },
                        grid: {
                            color: '#e0e0e0'
                        }
                    },
                    x: {
                        title: {
                            display: false
                        },
                        grid: {
                            color: '#e0e0e0'
                        }
                    }
                }
            }
        });
        
    } catch (error) {
        console.error(`[RJMCMC] Failed to create component chart ${canvasId}:`, error);
    }
}

// Control functions for MCMC state and burn-in toggle
function setMCMCRunningState(isRunning) {
    rjmcState.isRunning = isRunning;
    console.log(`[RJMCMC] MCMC running state set to: ${isRunning}`);
    
    // Update status text
    updatePlotStatus();
    
    // Update all plots with new running state
    updateComponentTracePlot();
    updateLogPosteriorPlot();
    // updateParameterTracePlots(); // Commented out for performance
}

function toggleBurninDisplay() {
    rjmcState.showBurnin = !rjmcState.showBurnin;
    console.log(`[RJMCMC] Burn-in display toggled to: ${rjmcState.showBurnin}`);
    
    // Update checkbox state
    const checkbox = document.getElementById('showBurninToggle');
    if (checkbox) {
        checkbox.checked = rjmcState.showBurnin;
    }
    
    // Update status text
    updatePlotStatus();
    
    // Update all plots with new burn-in setting
    updateComponentTracePlot();
    updateLogPosteriorPlot();
    // updateParameterTracePlots(); // Commented out for performance
}

function toggleDynamicParameterPlots() {
    rjmcState.dynamicParameterPlots = !rjmcState.dynamicParameterPlots;
    console.log(`[RJMCMC] Dynamic parameter plots toggled to: ${rjmcState.dynamicParameterPlots}`);
    
    // Update checkbox state
    const checkbox = document.getElementById('dynamicParameterPlotsToggle');
    if (checkbox) {
        checkbox.checked = rjmcState.dynamicParameterPlots;
    }
    
    // Show warning if enabling
    if (rjmcState.dynamicParameterPlots) {
        alert('WARNING: Dynamic parameter plot updates can significantly slow down the MCMC run. Only enable this if you need real-time parameter trace monitoring.');
    }
    
    // Update parameter plots immediately if enabling
    if (rjmcState.dynamicParameterPlots) {
        updateParameterTracePlots();
    }
}

// Manual function to test convergence diagnostics
function testConvergenceDiagnostics() {
    console.log('[RJMCMC] Manual convergence diagnostics test');
    updateConvergenceDiagnostics();
}

function updatePlotStatus() {
    const statusElement = document.getElementById('plotStatus');
    if (statusElement) {
        if (rjmcState.isRunning) {
            const maxPoints = rjmcState.maxDisplayPoints || 1000;
            statusElement.textContent = `Showing last ${maxPoints} points (MCMC running)`;
        } else {
            const burninText = rjmcState.showBurnin ? 'with burn-in' : 'without burn-in';
            statusElement.textContent = `Showing full trajectory (${burninText})`;
        }
    }
}

// Old mean trace plot function removed - now using faceted plots

// Old sigma trace plot function removed - now using faceted plots

// Get most common K value
function getMostCommonK() {
    let maxCount = 0;
    let mostCommonK = null;
    
    for (const k in rjmcState.kProbabilities) {
        if (rjmcState.kProbabilities[k] > maxCount) {
            maxCount = rjmcState.kProbabilities[k];
            mostCommonK = parseInt(k);
        }
    }
    
    return mostCommonK;
}

// Clean up old data to prevent memory buildup
function cleanupOldData() {
    // Clean up old posterior samples (keep only recent ones)
    const maxSamples = 10000; // Increased to 10000 for better posterior representation
    for (const k in rjmcState.kPosteriorSamples) {
        if (rjmcState.kPosteriorSamples[k].length > maxSamples) {
            rjmcState.kPosteriorSamples[k] = rjmcState.kPosteriorSamples[k].slice(-maxSamples);
        }
    }
    
    // Clean up chain histories to prevent memory bloat
    rjmcState.chains.forEach(chain => {
        const maxHistoryLength = 2000; // Keep only recent history
        
        if (chain.kHistory.length > maxHistoryLength) {
            chain.kHistory = chain.kHistory.slice(-maxHistoryLength);
        }
        if (chain.muHistory.length > maxHistoryLength) {
            chain.muHistory = chain.muHistory.slice(-maxHistoryLength);
        }
        if (chain.sigmaHistory.length > maxHistoryLength) {
            chain.sigmaHistory = chain.sigmaHistory.slice(-maxHistoryLength);
        }
        if (chain.logPostHistory.length > maxHistoryLength) {
            chain.logPostHistory = chain.logPostHistory.slice(-maxHistoryLength);
        }
    });
    
    // Reduce display points as run progresses to maintain performance
    if (rjmcState.iterations > 15000) {
        rjmcState.maxDisplayPoints = 500; // Show fewer points for very long runs
    } else if (rjmcState.iterations > 10000) {
        rjmcState.maxDisplayPoints = 750; // Show fewer points for long runs
    }
    
    // Force garbage collection hint (if available)
    if (typeof gc === 'function') {
        gc();
    }
    
    console.log('[RJMCMC] Cleaned up old data at iteration', rjmcState.iterations, 
                `(maxDisplayPoints: ${rjmcState.maxDisplayPoints})`);
}

// Update mixture fit
function updateMixtureFit() {
    console.log('[RJMCMC] updateMixtureFit called');
    
    if (!currentData) {
        console.log('[RJMCMC] No current data for mixture fit');
        return;
    }
    
    const canvas = document.getElementById('mixtureFitChart');
    if (!canvas) {
        console.error('[RJMCMC] mixtureFitChart canvas not found');
        return;
    }
    
    console.log('[RJMCMC] Canvas found:', canvas);
    console.log('[RJMCMC] kPosteriorSamples:', rjmcState.kPosteriorSamples);
    console.log('[RJMCMC] currentTab:', rjmcState.currentTab);
    
    const ctx = canvas.getContext('2d');
    
    if (mixtureFitChart) {
        console.log('[RJMCMC] Destroying existing mixture fit chart');
        mixtureFitChart.destroy();
        mixtureFitChart = null;
    }
    
    try {
        // Show raw data histogram
        const histogram = createHistogram(currentData.obs, 40);
        
        // Create fitted mixture density from posterior samples
        const k = parseInt(rjmcState.currentTab ? rjmcState.currentTab.replace('top', '') : '2');
        const samples = rjmcState.kPosteriorSamples[k] || [];
        
        console.log(`[RJMCMC] Creating mixture fit for K=${k} with ${samples.length} samples`);
        
        let fittedDensity;
        if (samples.length === 0) {
            // No samples yet, create zero density
            fittedDensity = histogram.x.map(() => 0);
        } else {
            fittedDensity = histogram.x.map(x => {
                let totalDensity = 0;
                
                // Average over all posterior samples for this K
                samples.forEach(sample => {
                    let sampleDensity = 0;
                    
                    // Order the means from lowest to highest for this sample
                    const muValues = sample.jump.mu.slice();
                    const sigmaValues = sample.jump.sigma.slice();
                    const pValues = sample.jump.p.slice();
                    
                    // Create indices sorted by mean values
                    const sortedIndices = muValues
                        .map((val, idx) => ({ val, idx }))
                        .sort((a, b) => a.val - b.val)
                        .map(item => item.idx);
                    
                    // Use ordered values
                    for (let j = 0; j < sortedIndices.length; j++) {
                        const originalIdx = sortedIndices[j];
                        sampleDensity += pValues[originalIdx] * normalDensity(parseFloat(x), muValues[originalIdx], sigmaValues[originalIdx]);
                    }
                    
                    totalDensity += sampleDensity;
                });
                
                // Average over all samples
                return (totalDensity / samples.length) * 100; // Scale for visualization
            });
        }
        
        // Create datasets array - only show fitted mixture, not raw data
        const datasets = [];
        
        // Only add fitted mixture if we have samples
        if (samples.length > 0) {
            datasets.push({
                label: `Fitted Mixture (${samples.length} samples)`,
                data: fittedDensity,
                type: 'line',
                borderColor: '#000000',
                backgroundColor: 'transparent',
                borderWidth: 3,
                fill: false,
                tension: 0.1
            });
        } else {
            // Show placeholder when no samples
            datasets.push({
                label: 'Fitted Mixture (No data yet)',
                data: histogram.x.map(() => 0),
                type: 'line',
                borderColor: '#cccccc',
                backgroundColor: 'transparent',
                borderWidth: 2,
                fill: false,
                borderDash: [5, 5]
            });
        }
        
        // Add vertical lines for component means if we have current jump data
        if (samples.length > 0 && rjmcState.currentJump && rjmcState.currentJump.mu) {
            const maxDensity = Math.max(...fittedDensity);
            const binWidth = parseFloat(histogram.x[1]) - parseFloat(histogram.x[0]);
            
            rjmcState.currentJump.mu.forEach((mean, index) => {
                // Find the bin index for the mean
                const meanIndex = histogram.x.findIndex(x => {
                    const binCenter = parseFloat(x);
                    return Math.abs(binCenter - mean) <= binWidth / 2;
                });
                
                if (meanIndex !== -1) {
                    // Create vertical line for mean
                    const meanLineData = new Array(histogram.x.length).fill(null);
                    meanLineData[meanIndex] = maxDensity * 0.9;
                
                    datasets.push({
                        label: `Component ${index + 1} Mean (${mean.toFixed(2)})`,
                        data: meanLineData,
                        type: 'line',
                        borderColor: `hsl(${index * 60}, 70%, 50%)`,
                        backgroundColor: `hsl(${index * 60}, 70%, 50%)`,
                        borderWidth: 4,
                        pointRadius: 6,
                        pointHoverRadius: 8,
                        fill: false,
                        tension: 0,
                        borderDash: [0, 0], // Solid line
                        showLine: true,
                        spanGaps: false
                    });
                }
            });
        }
        
        // Calculate appropriate y-axis maximum based on fitted density only
        const maxFittedDensity = samples.length > 0 ? Math.max(...fittedDensity) : 0;
        const yMax = maxFittedDensity > 0 ? maxFittedDensity * 1.2 : 10; // Default max if no data
        
        mixtureFitChart = new Chart(ctx, {
            type: 'line',
            data: {
                labels: histogram.x,
                datasets: datasets
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                height: 300, // Fixed height
                scales: {
                    y: {
                        beginAtZero: true,
                        max: yMax, // Set max based on fitted density only
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
                },
                plugins: {
                    title: {
                        display: true,
                        text: 'Mixture Model Fit'
                    },
                    legend: {
                        display: true,
                        position: 'top',
                        labels: {
                            boxWidth: 12,
                            padding: 8,
                            usePointStyle: true
                        }
                    }
                }
            }
        });
        
    } catch (error) {
        console.error('[RJMCMC] Failed to update mixture fit:', error);
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

// Test if Plotly is loaded
function testPlotlyLoaded() {
    console.log('[RJMCMC] Testing if Plotly is loaded...');
    console.log('[RJMCMC] typeof Plotly:', typeof Plotly);
    console.log('[RJMCMC] Plotly object:', Plotly);
    
    if (typeof Plotly === 'undefined') {
        console.error('[RJMCMC] Plotly is not loaded!');
        return false;
    } else {
        console.log('[RJMCMC] Plotly is loaded successfully');
        return true;
    }
}

// Simple test to create a basic plot
function createSimpleTestPlot() {
    console.log('[RJMCMC] Creating simple test plot...');
    
    const testData = [{
        x: [1, 2, 3, 4, 5],
        y: [2, 3, 2, 4, 3],
        type: 'scatter',
        mode: 'lines+markers',
        name: 'Test Line'
    }];
    
    const layout = {
        title: 'Simple Test Plot',
        xaxis: { title: 'X', autorange: true },
        yaxis: { title: 'Y', autorange: true, fixedrange: false }
    };
    
    try {
        Plotly.newPlot('rawDataChart', testData, layout);
        console.log('[RJMCMC] Simple test plot created successfully');
    } catch (error) {
        console.error('[RJMCMC] Simple test plot failed:', error);
    }
}

// Test function for raw data plot
function testRawDataPlot() {
    console.log('[RJMCMC] Testing raw data plot...');
    
    const testData = [{
        x: [1, 2, 2, 3, 3, 3, 4, 4, 5],
        type: 'histogram',
        nbinsx: 5,
        marker: { color: 'rgba(0,0,0,0.7)' },
        name: 'Test Data'
    }];
    
    const layout = {
        title: 'Test Raw Data',
        xaxis: { title: 'Value' },
        yaxis: { title: 'Frequency' }
    };
    
    try {
        Plotly.newPlot('rawDataChart', testData, layout);
        console.log('[RJMCMC] Test raw data plot created successfully');
    } catch (error) {
        console.error('[RJMCMC] Test raw data plot failed:', error);
    }
}

// Show raw data visualization
function showRawDataVisualization() {
    console.log('[RJMCMC] showRawDataVisualization called');
    console.log('[RJMCMC] currentData:', currentData);
    
    // Test if Plotly is loaded
    if (!testPlotlyLoaded()) {
        console.error('[RJMCMC] Cannot create plots - Plotly not loaded');
        return;
    }
    
    if (!currentData || currentData.obs.length === 0) {
        console.log('[RJMCMC] No currentData or empty obs, showing simple test plot instead');
        createSimpleTestPlot();
        return;
    }
    
    console.log('[RJMCMC] Data length:', currentData.obs.length);
    console.log('[RJMCMC] First 10 data points:', currentData.obs.slice(0, 10));
    
    // Calculate basic statistics
    const mean = currentData.obs.reduce((sum, val) => sum + val, 0) / currentData.obs.length;
    const variance = currentData.obs.reduce((sum, val) => sum + Math.pow(val - mean, 2), 0) / currentData.obs.length;
    const stdDev = Math.sqrt(variance);
    const min = Math.min(...currentData.obs);
    const max = Math.max(...currentData.obs);
    
    // Create Plotly histogram
    const plotData = [{
        x: currentData.obs,
        type: 'histogram',
        nbinsx: 40,
        marker: {
            color: 'rgba(128,128,128,0.4)',
            line: {
                color: 'rgba(128,128,128,0.8)',
                width: 1
            }
        },
        name: 'Data Distribution'
    }];
    
    console.log('[RJMCMC] Raw data plotData created:', plotData);
    
    const layout = {
        title: {
            text: 'Raw Data Distribution',
            font: { size: 14, family: 'Avenir, sans-serif' }
        },
        xaxis: {
            title: 'Value',
            showgrid: true,
            gridcolor: '#e0e0e0',
            font: { family: 'Avenir, sans-serif' },
            autorange: true,
            fixedrange: false
        },
        yaxis: {
            title: 'Frequency',
            showgrid: true,
            gridcolor: '#e0e0e0',
            font: { family: 'Avenir, sans-serif' },
            autorange: true,
            fixedrange: false
        },
        margin: { l: 60, r: 20, t: 40, b: 60 },
        plot_bgcolor: 'rgba(0,0,0,0)',
        paper_bgcolor: 'rgba(0,0,0,0)',
        font: { family: 'Avenir, sans-serif' },
        autosize: true,
        height: 300,
        width: '100%'
    };
    
    const config = {
        responsive: true,
        displayModeBar: false
    };
    
    const element = document.getElementById('rawDataChart');
    console.log('[RJMCMC] Raw data canvas element found:', element);
    
    if (!element) {
        console.error('[RJMCMC] Canvas element rawDataChart not found!');
        return;
    }
    
    try {
        Plotly.newPlot('rawDataChart', plotData, layout, config);
        console.log('[RJMCMC] Raw data plot created successfully');
    } catch (error) {
        console.error('[RJMCMC] Failed to create raw data Plotly chart:', error);
    }
}

// Switch between mixture fit tabs
function switchMixtureTab(tabName) {
    console.log(`[RJMCMC] Switching to tab: ${tabName}`);
    
    // Update active tab button
    document.querySelectorAll('.tab').forEach(btn => btn.classList.remove('active'));
    const tabButton = document.getElementById(`tab-${tabName}`);
    if (tabButton) {
        tabButton.classList.add('active');
    }
    
    // Update current tab
    rjmcState.currentTab = tabName;
    
    // Show/hide download button based on available data
    const downloadBtn = document.getElementById('downloadPosteriorsBtn');
    if (downloadBtn) {
        const k = parseInt(tabName.replace('top', ''));
        const hasData = rjmcState.kPosteriorSamples[k] && rjmcState.kPosteriorSamples[k].length > 0;
        downloadBtn.style.display = hasData ? 'block' : 'none';
    }
    
    // Update the mixture fit visualization
    updateMixtureFit();
}

// Update mixture fit tabs based on available K values
function updateMixtureTabs() {
    // Get available K values with samples, ordered by actual frequency (most common first)
    const availableKValues = Object.keys(rjmcState.kPosteriorSamples)
        .map(Number)
        .filter(kVal => rjmcState.kPosteriorSamples[kVal] && rjmcState.kPosteriorSamples[kVal].length > 0)
        .sort((a, b) => {
            // Sort by actual frequency from kProbabilities (descending), then by K value (ascending)
            const freqA = rjmcState.kProbabilities[a] || 0;
            const freqB = rjmcState.kProbabilities[b] || 0;
            if (freqA !== freqB) {
                return freqB - freqA; // Higher frequency first
            }
            return a - b; // Lower K first if same frequency
        });
    
    console.log(`[RJMCMC] Available K values: ${availableKValues.join(', ')}`);
    console.log(`[RJMCMC] K probabilities:`, rjmcState.kProbabilities);
    
    // Debug: Show the mode (most frequent K)
    if (availableKValues.length > 0) {
        const mode = availableKValues[0];
        const modeFreq = rjmcState.kProbabilities[mode] || 0;
        console.log(`[RJMCMC] Mode is K=${mode} with frequency ${modeFreq}`);
        
        // Debug: Show all frequencies for comparison
        console.log(`[RJMCMC] All frequencies:`, availableKValues.map(k => `K=${k}: ${rjmcState.kProbabilities[k] || 0}`).join(', '));
    }
    
    if (availableKValues.length === 0) {
        // No data yet, show default tabs
        return;
    }
    
    // Update the tab container
    const tabsHeader = document.querySelector('.tabs-header');
    if (!tabsHeader) {
        console.log('[RJMCMC] tabs-header not found, creating it');
        const tabsContainer = document.getElementById('mixtureTabs');
        if (tabsContainer) {
            const header = document.createElement('div');
            header.className = 'tabs-header';
            tabsContainer.appendChild(header);
            return updateMixtureTabs(); // Retry after creating header
        }
        return;
    }
    
    // Check if we need to recreate tabs (new K values or different order)
    const existingTabs = Array.from(tabsHeader.children).map(btn => {
        const k = parseInt(btn.id.replace('tab-top', ''));
        return k;
    }).sort((a, b) => a - b);
    
    const needsRecreation = existingTabs.length !== availableKValues.length || 
        !existingTabs.every((k, i) => k === availableKValues[i]);
    
    if (needsRecreation) {
        console.log(`[RJMCMC] Recreating tabs - existing: [${existingTabs.join(', ')}], new: [${availableKValues.join(', ')}]`);
        
        // Clear existing tabs
        tabsHeader.innerHTML = '';
        
        // Create new tabs based on available K values
        availableKValues.forEach((k, index) => {
            const samples = rjmcState.kPosteriorSamples[k];
            const frequency = rjmcState.kProbabilities[k] || 0;
            const totalSamples = rjmcState.iterations * rjmcState.numChains;
            const percentage = totalSamples > 0 ? ((frequency / totalSamples) * 100).toFixed(1) : '0.0';
            const isActive = index === 0; // First tab is active
            
            const button = document.createElement('button');
            button.className = `tab ${isActive ? 'active' : ''}`;
            button.id = `tab-top${k}`;
            button.textContent = `K=${k} (${percentage}%)`;
            button.onclick = () => switchMixtureTab(`top${k}`);
            
            console.log(`[RJMCMC] Creating tab for K=${k} with ${frequency}/${totalSamples} samples (${percentage}%)`);
            console.log(`[RJMCMC] Tab text content: "${button.textContent}"`);
            tabsHeader.appendChild(button);
        });
        
        // Update current tab to the first available one
        if (availableKValues.length > 0) {
            const firstK = availableKValues[0];
            rjmcState.currentTab = `top${firstK}`;
            console.log(`[RJMCMC] Set current tab to K=${firstK}`);
            
            // Verify the first tab is actually the mode
            const firstFreq = rjmcState.kProbabilities[firstK] || 0;
            console.log(`[RJMCMC] First tab K=${firstK} has frequency ${firstFreq}`);
            
            // Show download button if data is available
            const downloadBtn = document.getElementById('downloadPosteriorsBtn');
            if (downloadBtn) {
                const hasData = rjmcState.kPosteriorSamples[firstK] && rjmcState.kPosteriorSamples[firstK].length > 0;
                downloadBtn.style.display = hasData ? 'block' : 'none';
            }
        }
    } else {
        // Just update the sample counts in existing tabs
        availableKValues.forEach((k, index) => {
            const samples = rjmcState.kPosteriorSamples[k];
            const tabButton = document.getElementById(`tab-top${k}`);
            if (tabButton) {
                tabButton.textContent = `K=${k} (${samples.length})`;
                console.log(`[RJMCMC] Updated tab K=${k} to show ${samples.length} samples`);
            }
        });
    }
}

// Update tab sample counts without recreating tabs
function updateTabSampleCounts() {
    const tabsHeader = document.querySelector('.tabs-header');
    if (!tabsHeader) return;
    
    // Update percentages for existing tabs
    Array.from(tabsHeader.children).forEach(btn => {
        const k = parseInt(btn.id.replace('tab-top', ''));
        const frequency = rjmcState.kProbabilities[k] || 0;
        const totalSamples = rjmcState.iterations * rjmcState.numChains;
        const percentage = totalSamples > 0 ? ((frequency / totalSamples) * 100).toFixed(1) : '0.0';
        
        if (frequency > 0) {
            btn.textContent = `K=${k} (${percentage}%)`;
            console.log(`[RJMCMC] Updated tab K=${k} text to: "${btn.textContent}"`);
        }
    });
}

// Force update tabs to show the mode first (for final results)
function forceUpdateTabsToShowMode() {
    console.log('[RJMCMC] Force updating tabs to show mode first');
    console.log('[RJMCMC] Final K probabilities:', rjmcState.kProbabilities);
    
    // Find the mode (most frequent K)
    let mode = null;
    let maxFreq = 0;
    console.log('[RJMCMC] Searching for mode in kProbabilities:', rjmcState.kProbabilities);
    
    for (const k in rjmcState.kProbabilities) {
        const freq = rjmcState.kProbabilities[k];
        console.log(`[RJMCMC] K=${k}: frequency=${freq}`);
        if (freq > maxFreq) {
            maxFreq = freq;
            mode = parseInt(k);
            console.log(`[RJMCMC] New mode candidate: K=${mode} with frequency=${maxFreq}`);
        }
    }
    
    if (mode !== null) {
        console.log(`[RJMCMC] Mode is K=${mode} with frequency ${maxFreq}`);
        
        // Force recreate tabs with mode first by manually ordering
        const tabsHeader = document.querySelector('.tabs-header');
        if (tabsHeader) {
            // Clear existing tabs completely
            tabsHeader.innerHTML = '';
            console.log('[RJMCMC] Cleared all existing tabs');
            
            // Get all K values with samples, but put mode first
            const allKValues = Object.keys(rjmcState.kPosteriorSamples)
                .map(Number)
                .filter(kVal => rjmcState.kPosteriorSamples[kVal] && rjmcState.kPosteriorSamples[kVal].length > 0);
            
            // Put mode first, then sort the rest by frequency
            const otherKValues = allKValues.filter(k => k !== mode).sort((a, b) => {
                const freqA = rjmcState.kProbabilities[a] || 0;
                const freqB = rjmcState.kProbabilities[b] || 0;
                return freqB - freqA;
            });
            
            const orderedKValues = [mode, ...otherKValues];
            console.log(`[RJMCMC] Ordered K values: ${orderedKValues.join(', ')}`);
            
            // Create tabs in the correct order
            orderedKValues.forEach((k, index) => {
                const samples = rjmcState.kPosteriorSamples[k];
                const frequency = rjmcState.kProbabilities[k] || 0;
                const totalSamples = rjmcState.iterations * rjmcState.numChains;
                const percentage = totalSamples > 0 ? ((frequency / totalSamples) * 100).toFixed(1) : '0.0';
                const isActive = index === 0; // First tab (mode) is active
                
                const button = document.createElement('button');
                button.className = `tab ${isActive ? 'active' : ''}`;
                button.id = `tab-top${k}`;
                button.textContent = `K=${k} (${percentage}%)`;
                button.onclick = () => switchMixtureTab(`top${k}`);
                
                console.log(`[RJMCMC] Creating tab for K=${k} with ${frequency}/${totalSamples} samples (${percentage}%)`);
                console.log(`[RJMCMC] Tab text content: "${button.textContent}"`);
                tabsHeader.appendChild(button);
            });
            
            // Set current tab to mode
            rjmcState.currentTab = `top${mode}`;
            console.log(`[RJMCMC] Set current tab to mode K=${mode}`);
            
            // Show download button if data is available
            const downloadBtn = document.getElementById('downloadPosteriorsBtn');
            if (downloadBtn) {
                const hasData = rjmcState.kPosteriorSamples[mode] && rjmcState.kPosteriorSamples[mode].length > 0;
                downloadBtn.style.display = hasData ? 'block' : 'none';
            }
        }
    } else {
        console.log('[RJMCMC] No mode found, using regular tab update');
        updateMixtureTabs();
    }
}

// Set up automatic tab switching to show available data
function setupAutoTabSwitching() {
    // Check every 100 iterations if we should update tabs and switch to one with data
    const checkInterval = setInterval(() => {
        if (rjmcState.iterations > 0) {
            // Update tabs based on available data
            updateMixtureTabs();
            
            // Update the mixture fit with current tab
            updateMixtureFit();
        }
    }, 1000); // Check every second
    
    // Store the interval ID so we can clear it later if needed
    rjmcState.autoTabCheckInterval = checkInterval;
}

// Clean up auto tab switching
function cleanupAutoTabSwitching() {
    if (rjmcState.autoTabCheckInterval) {
        clearInterval(rjmcState.autoTabCheckInterval);
        rjmcState.autoTabCheckInterval = null;
    }
}




// Initialize file upload functionality
function initializeFileUpload() {
    const fileUpload = document.getElementById('fileUpload');
    const fileInput = document.getElementById('fileInput');
    const fileInfo = document.getElementById('fileInfo');

    if (!fileUpload || !fileInput) {
        console.warn('[RJMCMC] File upload elements not found');
        return;
    }

    fileUpload.addEventListener('click', () => fileInput.click());
    fileUpload.addEventListener('dragover', handleDragOver);
    fileUpload.addEventListener('dragleave', handleDragLeave);
    fileUpload.addEventListener('drop', handleDrop);
    fileInput.addEventListener('change', handleFileSelect);
    
    console.log('[RJMCMC] File upload functionality initialized');
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

    const reader = new FileReader();
    reader.onload = function(e) {
        try {
            const csvText = e.target.result;
            const lines = csvText.split('\n');
            
            // Extract numeric data from first column, skip empty lines
            const data = lines
                .map(line => line.trim())
                .filter(line => line.length > 0)
                .map(line => {
                    // Split by comma and take first column
                    const columns = line.split(',');
                    return parseFloat(columns[0]);
                })
                .filter(val => !isNaN(val));

            if (data.length === 0) {
                showError('No valid numeric data found in CSV file.');
                return;
            }

            // Update current data
            currentData = {
                obs: data,
                n: data.length
            };
            
            // Update mean range based on data
            updateMeanRangeFromData(data);
            
            document.getElementById('fileInfo').innerHTML = 
                `Loaded ${data.length} data points from ${file.name}`;
            
            showSuccess(`Successfully loaded ${data.length} data points.`);
            
            // Show raw data visualization
            showRawDataVisualization();
            
            console.log('[RJMCMC] CSV file loaded:', data.length, 'data points');
        } catch (error) {
            showError('Error reading file: ' + error.message);
        }
    };
    
    reader.readAsText(file);
}

function showError(message) {
    console.error('[RJMCMC] Error:', message);
    // You could add a visual error display here
    alert('Error: ' + message);
}

function showSuccess(message) {
    console.log('[RJMCMC] Success:', message);
    // You could add a visual success display here
}

// Update mean range inputs based on data
function updateMeanRangeFromData(data) {
    const minVal = Math.min(...data);
    const maxVal = Math.max(...data);
    const range = maxVal - minVal;
    
    // Add some padding (20% of range on each side)
    const padding = range * 0.2;
    const minMean = minVal - padding;
    const maxMean = maxVal + padding;
    
    // Update the mean range inputs
    const minMeanInput = document.getElementById('minMean');
    const maxMeanInput = document.getElementById('maxMean');
    
    if (minMeanInput) {
        minMeanInput.value = minMean.toFixed(2);
        minMeanInput.min = minMean - padding; // Allow even lower values
        minMeanInput.max = maxMean + padding; // Allow even higher values
    }
    
    if (maxMeanInput) {
        maxMeanInput.value = maxMean.toFixed(2);
        maxMeanInput.min = minMean - padding; // Allow even lower values
        maxMeanInput.max = maxMean + padding; // Allow even higher values
    }
    
    console.log(`[RJMCMC] Updated mean range: [${minMean.toFixed(2)}, ${maxMean.toFixed(2)}] based on data range [${minVal.toFixed(2)}, ${maxVal.toFixed(2)}]`);
}

// Load sample data - nice 3-component example for vignette
function loadSampleData() {
    // Generate a nice 3-component mixture example
    const nSamples = 200;
    const trueK = 3;
    
    // True parameters for a nice 3-component mixture
    const trueMeans = [-2.0, 1.0, 4.0];
    const trueSigmas = [0.8, 0.6, 1.0];
    const trueMixingProps = [0.3, 0.4, 0.3];
    
    // Generate data
    const data = [];
    for (let i = 0; i < nSamples; i++) {
        // Choose component
        const u = Math.random();
        let component = 0;
        let cumProb = trueMixingProps[0];
        
        for (let k = 0; k < trueK - 1; k++) {
            if (u < cumProb) {
                component = k;
                break;
            }
            cumProb += trueMixingProps[k + 1];
            component = k + 1;
        }
        
        // Generate observation from chosen component using Box-Muller transform
        const mean = trueMeans[component];
        const sigma = trueSigmas[component];
        
        // Box-Muller transform for normal distribution
        const u1 = Math.random();
        const u2 = Math.random();
        const z0 = Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
        const obs = mean + sigma * z0;
        
        data.push(obs);
    }
    
    // Update current data
    currentData = {
        obs: data,
        n: data.length
    };
    
    // Update mean range based on data
    updateMeanRangeFromData(data);
    
    document.getElementById('fileInfo').innerHTML = 
        `Loaded ${data.length} data points from 3-component mixture example`;
    
    showSuccess(`Successfully loaded ${data.length} data points from 3-component mixture example.`);
    
    // Show raw data visualization
    showRawDataVisualization();
    
    console.log('[RJMCMC] Sample data loaded:', data.length, 'data points from 3-component mixture example');
    console.log('[RJMCMC] True parameters - Means:', trueMeans, 'Sigmas:', trueSigmas, 'Mixing props:', trueMixingProps);
}

// Stop RJMCMC analysis
function stopRJMCMCAnalysis() {
    if (rjmcState.isAnalysisRunning) {
        rjmcState.isAnalysisRunning = false;
        console.log('[RJMCMC] Analysis stopped by user');
        
        // Update UI
        document.getElementById('startBtn').disabled = false;
        document.getElementById('stopBtn').disabled = true;
        
        const status = document.getElementById('status');
        status.className = 'status error';
        status.textContent = 'Analysis stopped by user';
        status.style.display = 'block';
        
        // Hide progress bar
        hideProgressBar();
    }
}

// Initialize the application
async function initSimpleApp() {
    try {
        console.log('[RJMCMC] Initializing simplified RJMCMC application...');
        
        // Initialize file upload functionality
        initializeFileUpload();
        
        // Generate initial sample data
        currentData = generateSampleData();
        console.log('[RJMCMC] Generated sample data:', currentData.obs.length, 'points');
        
        // Show raw data visualization
        showRawDataVisualization();
        
        // Initialize trace plots with empty data
        updateComponentTracePlot();
        updateLogPosteriorPlot();
        
        // Component distribution chart will be shown at the end of the run
        
        // Initialize mixture fit with raw data
        updateMixtureFit();
        
        // Initialize component distribution
        updateComponentDistribution();
        
        // Show posterior distribution container
        const posteriorContainer = document.getElementById('posteriorDistContainer');
        if (posteriorContainer) {
            posteriorContainer.style.display = 'block';
            console.log('[RJMCMC] Posterior distribution container made visible during init');
        }
        
        
        // Show results container
        const resultsContainer = document.getElementById('resultsContainer');
        const welcomeMessage = document.getElementById('welcomeMessage');
        
        if (resultsContainer) {
            resultsContainer.style.display = 'block';
        }
        if (welcomeMessage) {
            welcomeMessage.style.display = 'none';
        }
        
        // Enable analysis button
        const runAnalysisBtn = document.getElementById('runAnalysis');
        if (runAnalysisBtn) {
            runAnalysisBtn.disabled = false;
            runAnalysisBtn.textContent = 'Start RJMCMC Analysis';
        }
        
        console.log('[RJMCMC] Application initialized successfully');
        
    } catch (error) {
        console.error('[RJMCMC] Failed to initialize application:', error);
        const errorMessage = document.getElementById('error-message');
        const errorText = document.getElementById('error-text');
        if (errorMessage) {
            errorMessage.style.display = 'block';
        }
        if (errorText) {
            errorText.textContent = `Initialization Error: ${error.message}`;
        }
    }
}

// Start RJMCMC analysis
async function startRJMCMCAnalysis() {
    try {
        console.log('[RJMCMC] Starting RJMCMC analysis...');
        
        // Check if data is available
        if (!currentData || !currentData.obs || currentData.obs.length === 0) {
            throw new Error('No data available. Please load demo data first.');
        }
        
        // Get settings from UI
        const iterationsInput = document.getElementById('iterations');
        const burninInput = document.getElementById('burnin');
        
        const iterations = iterationsInput ? parseInt(iterationsInput.value) || 20000 : 20000;
        const burnin = burninInput ? parseInt(burninInput.value) || 10000 : 10000;
        
        // Get prior configuration from UI
        const poissonLambdaInput = document.getElementById('poissonLambda');
        const regularizationStrengthInput = document.getElementById('regularizationStrength');
        
        const poissonLambda = poissonLambdaInput ? parseFloat(poissonLambdaInput.value) || 2.0 : 2.0;
        const regularizationStrength = regularizationStrengthInput ? parseFloat(regularizationStrengthInput.value) || 0.1 : 0.1;
        
        // Update prior configuration
        PRIOR_CONFIG.poissonLambda = poissonLambda;
        PRIOR_CONFIG.regularizationStrength = regularizationStrength;
        
        console.log('[RJMCMC] Settings:', { iterations, burnin, poissonLambda, regularizationStrength });
        
        const settings = {
            iterations: iterations,
            burnin: burnin,
            updateInterval: 50
        };
        
        // Set analysis running state
        rjmcState.isAnalysisRunning = true;
        
        // Disable button during analysis
        const startBtn = document.getElementById('startBtn');
        const stopBtn = document.getElementById('stopBtn');
        if (startBtn) {
            startBtn.disabled = true;
            startBtn.textContent = 'Running Analysis...';
        }
        if (stopBtn) {
            stopBtn.disabled = false;
        }
        
        // Update status
        const status = document.getElementById('status');
        if (status) {
            status.className = 'status running';
            status.textContent = 'Starting RJMCMC analysis...';
            status.style.display = 'block';
        }
        
        // Run analysis
        await runRJMCMCAnalysis(settings);
        
        // Analysis completed
        rjmcState.isAnalysisRunning = false;
        
        // Re-enable button
        if (startBtn) {
            startBtn.disabled = false;
            startBtn.textContent = 'Start Analysis';
        }
        if (stopBtn) {
            stopBtn.disabled = true;
        }
        
        // Update status
        if (status) {
            status.className = 'status completed';
            status.textContent = 'Analysis completed successfully!';
        }
        
    } catch (error) {
        console.error('[RJMCMC] RJMCMC analysis failed:', error);
        
        // Reset analysis state
        rjmcState.isAnalysisRunning = false;
        
        // Re-enable button
        const startBtn = document.getElementById('startBtn');
        const stopBtn = document.getElementById('stopBtn');
        if (startBtn) {
            startBtn.disabled = false;
            startBtn.textContent = 'Start Analysis';
        }
        if (stopBtn) {
            stopBtn.disabled = true;
        }
        
        // Show error
        const status = document.getElementById('status');
        if (status) {
            status.className = 'status error';
            status.textContent = `Analysis Error: ${error.message}`;
            status.style.display = 'block';
        }
    }
}

// Load demo data
function loadDemoData() {
    try {
        console.log('[RJMCMC] Loading demo data...');
        currentData = generateSampleData();
        console.log('[RJMCMC] Demo data loaded:', currentData.obs.length, 'points');
        
        showRawDataVisualization();
        
        const resultsContainer = document.getElementById('resultsContainer');
        const welcomeMessage = document.getElementById('welcomeMessage');
        const runAnalysisBtn = document.getElementById('runAnalysis');
        
        if (resultsContainer) {
            resultsContainer.style.display = 'block';
        }
        if (welcomeMessage) {
            welcomeMessage.style.display = 'none';
        }
        if (runAnalysisBtn) {
            runAnalysisBtn.disabled = false;
        }
        
        console.log('[RJMCMC] Demo data setup complete');
    } catch (error) {
        console.error('[RJMCMC] Failed to load demo data:', error);
    }
}


// Download posteriors function
function downloadPosteriors() {
    const k = parseInt(rjmcState.currentTab ? rjmcState.currentTab.replace('top', '') : '2');
    const samples = rjmcState.kPosteriorSamples[k] || [];
    
    if (samples.length === 0) {
        alert('No posterior samples available for download');
        return;
    }
    
    // Convert samples to CSV format
    let csvContent = 'Sample,Component,Mean,Sigma,Proportion\n';
    
    samples.forEach((sample, sampleIndex) => {
        const muValues = sample.jump.mu.slice();
        const sigmaValues = sample.jump.sigma.slice();
        const pValues = sample.jump.p.slice();
        
        // Order by mean values
        const sortedIndices = muValues
            .map((val, idx) => ({ val, idx }))
            .sort((a, b) => a.val - b.val)
            .map(item => item.idx);
        
        sortedIndices.forEach((originalIdx, compIndex) => {
            csvContent += `${sampleIndex + 1},${compIndex + 1},${muValues[originalIdx]},${sigmaValues[originalIdx]},${pValues[originalIdx]}\n`;
        });
    });
    
    // Create and download file
    const blob = new Blob([csvContent], { type: 'text/csv' });
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `posterior_samples_K${k}.csv`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    window.URL.revokeObjectURL(url);
}

// Progress bar functions
function showProgressBar() {
    const progressContainer = document.getElementById('progressContainer');
    if (progressContainer) {
        progressContainer.style.display = 'block';
    }
}

function hideProgressBar() {
    const progressContainer = document.getElementById('progressContainer');
    if (progressContainer) {
        progressContainer.style.display = 'none';
    }
}

function updateProgressBar(currentIteration, totalIterations) {
    const progressBar = document.getElementById('progressBar');
    const progressPercent = document.getElementById('progressPercent');
    const currentIterationSpan = document.getElementById('currentIteration');
    const totalIterationsSpan = document.getElementById('totalIterations');
    
    if (progressBar && progressPercent && currentIterationSpan && totalIterationsSpan) {
        const percentage = Math.min(100, Math.max(0, (currentIteration / totalIterations) * 100));
        
        progressBar.style.width = percentage + '%';
        progressPercent.textContent = Math.round(percentage) + '%';
        currentIterationSpan.textContent = currentIteration.toLocaleString();
        totalIterationsSpan.textContent = totalIterations.toLocaleString();
    }
}

// Test function for posterior distribution
function testPosteriorDistribution() {
    console.log('[RJMCMC] Testing posterior distribution...');
    
    // Create test data for all chains
    console.log('[RJMCMC] Creating test chains data');
    rjmcState.chains = [
        {
            kHistory: [2, 2, 2, 3, 3, 2, 3, 3, 4, 3, 2, 3, 4, 4, 3, 2, 3, 2, 3, 3, 2, 2, 3, 4, 3, 2, 3, 3, 4, 3]
        },
        {
            kHistory: [2, 3, 3, 3, 2, 3, 4, 3, 3, 2, 3, 4, 3, 2, 3, 3, 4, 3, 2, 3, 3, 2, 3, 4, 3, 2, 3, 3, 4, 3]
        },
        {
            kHistory: [2, 2, 3, 2, 3, 3, 4, 3, 2, 3, 4, 4, 3, 2, 3, 3, 2, 3, 4, 3, 2, 3, 3, 4, 3, 2, 3, 3, 4, 3]
        },
        {
            kHistory: [3, 3, 2, 3, 4, 3, 2, 3, 3, 4, 3, 2, 3, 3, 2, 3, 4, 3, 3, 2, 3, 4, 3, 2, 3, 3, 4, 3, 2, 3]
        }
    ];
    
    console.log('[RJMCMC] Test chains created:', rjmcState.chains.length);
    updateComponentDistribution();
}

// Test function for mixture fit
function testMixtureFit() {
    console.log('[RJMCMC] Testing mixture fit...');
    
    // Create test data
    if (!currentData) {
        console.log('[RJMCMC] Creating test data for mixture fit');
        currentData = {
            obs: Array.from({length: 100}, () => Math.random() * 10 - 5)
        };
    }
    
    // Create test posterior samples
    if (!rjmcState.kPosteriorSamples || Object.keys(rjmcState.kPosteriorSamples).length === 0) {
        console.log('[RJMCMC] Creating test posterior samples');
        rjmcState.kPosteriorSamples = {
            '2': Array.from({length: 50}, () => ({
                jump: {
                    p: [0.5, 0.5],
                    mu: [1, 4],
                    sigma: [1, 1]
                },
                logPosterior: -100
            })),
            '3': Array.from({length: 30}, () => ({
                jump: {
                    p: [0.33, 0.33, 0.34],
                    mu: [0, 2, 5],
                    sigma: [1, 1, 1]
                },
                logPosterior: -95
            }))
        };
        rjmcState.currentTab = 'top2';
    }
    
    console.log('[RJMCMC] Test data created, calling updateMixtureFit');
    updateMixtureFit();
}

// Make functions globally available
window.initSimpleApp = initSimpleApp;
window.startRJMCMCAnalysis = startRJMCMCAnalysis;
window.showRawDataVisualization = showRawDataVisualization;
window.testPosteriorDistribution = testPosteriorDistribution;
window.testMixtureFit = testMixtureFit;
window.testPlotly = createSimpleTestPlot;
window.testComponentPlot = testComponentTracePlot;
window.loadDemoData = loadDemoData;
window.downloadPosteriors = downloadPosteriors;

// Add a simple test function
window.testButton = function() {
    console.log('[RJMCMC] Test button clicked!');
    alert('Button click is working!');
};

// Initialize when DOM is loaded and Plotly is available
function waitForPlotlyAndInit() {
    if (typeof Plotly !== 'undefined') {
        console.log('[RJMCMC] Plotly is loaded, initializing app...');
    initSimpleApp();
    } else {
        console.log('[RJMCMC] Waiting for Plotly to load...');
        setTimeout(waitForPlotlyAndInit, 100);
    }
}

document.addEventListener('DOMContentLoaded', function() {
    console.log('[RJMCMC] DOM loaded, waiting for Plotly...');
    waitForPlotlyAndInit();
});

// Also try to initialize immediately in case DOM is already loaded
if (document.readyState === 'loading') {
    console.log('[RJMCMC] DOM still loading...');
} else {
    console.log('[RJMCMC] DOM already loaded, waiting for Plotly...');
    waitForPlotlyAndInit();
}
