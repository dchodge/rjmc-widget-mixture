# Mathematical Foundations of Reversible Jump Markov Chain Monte Carlo for Mixture Distributions

## Introduction

This document provides a comprehensive mathematical treatment of the Reversible Jump Markov Chain Monte Carlo (RJMCMC) algorithm for fitting mixture distributions, as implemented in the interactive widget. The problem of determining the number of components in a mixture model is a classic example of model selection in Bayesian statistics, where the parameter space has varying dimensionality.

### Problem Statement

Given a dataset **y** = {y₁, y₂, ..., yₙ}, we wish to fit a finite mixture of normal distributions:

```
f(y|θ, K) = Σ(k=1 to K) πₖ N(y|μₖ, σₖ²)
```

where:
- K is the number of components (unknown)
- πₖ are the mixing proportions with Σ(k=1 to K) πₖ = 1
- μₖ and σₖ² are the mean and variance of component k
- θ = {π₁, ..., πₖ, μ₁, ..., μₖ, σ₁², ..., σₖ²}

The challenge is to simultaneously estimate both the number of components K and the parameters θ, where the dimensionality of θ depends on K.

## Mathematical Framework

### Likelihood Function

For a dataset **y** = {y₁, y₂, ..., yₙ}, the likelihood function is:

```
L(θ, K|y) = ∏(i=1 to n) Σ(k=1 to K) πₖ N(yᵢ|μₖ, σₖ²)
```

Taking the logarithm:

```
ℓ(θ, K|y) = Σ(i=1 to n) log(Σ(k=1 to K) πₖ N(yᵢ|μₖ, σₖ²))
```

where the normal density is:

```
N(y|μ, σ²) = (1/√(2πσ²)) exp(-(y-μ)²/(2σ²))
```

### Prior Distributions

We specify the following prior distributions:

#### Prior for the Number of Components

```
p(K) = 1/(K_max - K_min + 1), K ∈ {K_min, K_min+1, ..., K_max}
```

In our implementation: K_min = 1, K_max = 8.

#### Prior for Mixing Proportions

For each K, the mixing proportions follow a Dirichlet distribution:

```
π = (π₁, ..., πₖ) ~ Dirichlet(α₁, ..., αₖ)
```

With uniform priors: αₖ = 1 for all k.

#### Prior for Component Means

```
μₖ ~ Uniform(a, b), k = 1, ..., K
```

In our implementation: a = 0, b = 20.

#### Prior for Component Variances

```
σₖ² ~ Uniform(c, d), k = 1, ..., K
```

In our implementation: c = 0.3, d = 3.

### Posterior Distribution

The joint posterior distribution is:

```
p(K, θ|y) ∝ L(θ, K|y) · p(θ|K) · p(K)
```

The log-posterior is:

```
log p(K, θ|y) = ℓ(θ, K|y) + log p(θ|K) + log p(K) + constant
```

## Reversible Jump Markov Chain Monte Carlo

RJMCMC allows sampling from distributions with varying dimensionality by constructing a Markov chain that can jump between different model spaces.

### Proposal Types

At each iteration, we propose one of three moves:

1. **Birth move**: Increase K by 1
2. **Death move**: Decrease K by 1  
3. **Jump move**: Update parameters within the current model

### Birth Move

#### Proposal Mechanism

When proposing a birth move from model K to model K+1:

1. Generate new mixing proportion: π_{K+1} ~ Beta(1, K)
2. Renormalize existing proportions: πₖ' = πₖ(1 - π_{K+1}) for k = 1, ..., K
3. Generate new component parameters:
   - μ_{K+1} ~ Uniform(a, b)
   - σ_{K+1}² ~ Uniform(c, d)

#### Acceptance Probability

The acceptance probability for a birth move is:

```
α_birth = min(1, (p(K+1, θ'|y)/p(K, θ|y)) · (q_death(θ|θ')/q_birth(θ'|θ)) · (p_death/p_birth))
```

where:
- q_birth(θ'|θ) is the proposal density for birth
- q_death(θ|θ') is the proposal density for death
- p_birth and p_death are the probabilities of proposing birth and death moves

#### Jacobian

The Jacobian for the birth move transformation is:

```
J_birth = ∂(π₁', ..., πₖ', π_{K+1}, μ_{K+1}, σ_{K+1}²)/∂(π₁, ..., πₖ, u₁, u₂, u₃)
```

where u₁, u₂, u₃ are the random variables used to generate the new component.

### Death Move

#### Proposal Mechanism

When proposing a death move from model K to model K-1:

1. Select component j to remove uniformly: j ~ Uniform{1, ..., K}
2. Redistribute the removed proportion: πₖ' = πₖ/(1 - πⱼ) for k ≠ j
3. Remove μⱼ and σⱼ²

#### Acceptance Probability

The acceptance probability for a death move is:

```
α_death = min(1, (p(K-1, θ'|y)/p(K, θ|y)) · (q_birth(θ|θ')/q_death(θ'|θ)) · (p_birth/p_death))
```

### Jump Move

For jump moves within the same model, we use standard Metropolis-Hastings updates:

1. Update mixing proportions: π' ~ Dirichlet(α π)
2. Update component means: μₖ' ~ N(μₖ, τ_μ²)
3. Update component variances: log σₖ' ~ N(log σₖ, τ_σ²)

## Implementation Details

### Log-Posterior Calculation

In the widget implementation, the log-posterior is calculated as:

```
log p(K, θ|y) = Σ(i=1 to n) log(Σ(k=1 to K) πₖ N(yᵢ|μₖ, σₖ²)) + log p(θ|K) + log p(K)
```

where:

```
log p(θ|K) = Σ(k=1 to K) [log p(μₖ) + log p(σₖ²)] + log p(π)
           = Σ(k=1 to K) [log(1/(b-a)) + log(1/(d-c))] + log(Γ(Σ(k=1 to K) αₖ)/∏(k=1 to K) Γ(αₖ) ∏(k=1 to K) πₖ^(αₖ-1))
```

### Multiple Chains

The widget runs M = 4 parallel chains to assess convergence and provide robust estimates. Each chain m maintains its own state:

```
θ^(m) = {π₁^(m), ..., π_{K^(m)}^(m), μ₁^(m), ..., μ_{K^(m)}^(m), σ₁^(m), ..., σ_{K^(m)}^(m)}
```

### Posterior Summaries

#### Component Probabilities

The posterior probability of K components is estimated as:

```
p̂(K|y) = (1/(M·T)) Σ(m=1 to M) Σ(t=1 to T) I(K^(m,t) = K)
```

where T is the number of iterations and I(·) is the indicator function.

#### Parameter Estimates

For a given K, the posterior mean of parameter θₖ is:

```
θ̂ₖ = (1/(M·T_K)) Σ(m=1 to M) Σ(t: K^(m,t) = K) θₖ^(m,t)
```

where T_K is the number of iterations where K^(m,t) = K.

## Convergence Diagnostics

### Gelman-Rubin Statistic

For scalar parameters, the Gelman-Rubin statistic is:

```
R̂ = √(V̂/W)
```

where:
- W = (1/M) Σ(m=1 to M) s_m² (within-chain variance)
- B = (T/(M-1)) Σ(m=1 to M) (θ̄_m - θ̄)² (between-chain variance)
- V̂ = ((T-1)/T) W + (1/T) B (pooled variance)

### Effective Sample Size

The effective sample size is:

```
ESS = (MT)/(1 + 2Σ(t=1 to ∞) ρ_t)
```

where ρ_t is the autocorrelation at lag t.

## Widget Implementation

### Data Generation

The widget generates synthetic data from a 3-component mixture:

```
yᵢ ~ 0.3 · N(2, 1) + 0.7 · N(6, 0.25) + 1.0 · N(12, 1)
   = Σ(k=1 to 3) πₖ N(μₖ, σₖ²)
```

where π = (0.3, 0.7, 1.0), μ = (2, 6, 12), and σ² = (1, 0.25, 1).

### Visualization

The widget provides real-time visualization of:

1. **Trace plots**: Component count and log-posterior over iterations
2. **Component distribution**: Posterior probabilities for each K
3. **Mixture fit**: Data histogram with fitted density and component means
4. **Multiple chains**: All 4 chains displayed simultaneously

### Tabbed Interface

The mixture fit visualization includes tabs for:
- **Current K**: Real-time fit from the current state
- **Top K=X**: Best fit found for each number of components

Each tab shows the fitted mixture density with vertical lines indicating component means.

## Conclusion

The RJMCMC algorithm provides a principled approach to Bayesian model selection for mixture distributions. The mathematical framework ensures detailed balance and proper exploration of the model space, while the multiple-chain implementation provides robust convergence diagnostics. The interactive widget demonstrates these concepts through real-time visualization of the sampling process.

The key mathematical components include:
- Proper prior specification across varying dimensions
- Careful proposal mechanisms for birth/death moves
- Correct Jacobian calculations for dimension-changing moves
- Multiple-chain implementation for convergence assessment
- Real-time posterior summarization and visualization

This implementation serves as both an educational tool and a practical demonstration of advanced Bayesian computational methods.

## Mathematical Notation Summary

| Symbol | Meaning |
|--------|---------|
| K | Number of components |
| πₖ | Mixing proportion for component k |
| μₖ | Mean of component k |
| σₖ² | Variance of component k |
| θ | Parameter vector |
| y | Observed data |
| L(·) | Likelihood function |
| p(·) | Prior/posterior probability |
| α | Acceptance probability |
| q(·) | Proposal density |
| J | Jacobian determinant |
| M | Number of chains |
| T | Number of iterations |
| ESS | Effective sample size |
| R̂ | Gelman-Rubin statistic |
