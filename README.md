# RJMCMC Mixture Model Widget

A web-based interactive widget for dynamically determining the number of components in a mixture model using Reversible Jump Markov Chain Monte Carlo (RJMCMC).

## 🌐 Live Demo

The widget is automatically deployed to GitHub Pages: **[View Live Demo](https://dchodge.github.io/rjmc_widget/)**

## 🚀 Features

- **Dynamic Component Detection**: Automatically determines the optimal number of mixture components without guessing
- **Interactive Visualization**: Real-time plots showing MCMC sampling progress
- **CSV Data Upload**: Upload your own data or use provided examples
- **Convergence Diagnostics**: ESS, Gelman-Rubin statistic, and acceptance rates
- **Posterior Analysis**: Detailed posterior distributions for mixture parameters

## 📊 How It Works

The application uses **Reversible Jump MCMC** to explore different numbers of mixture components (K) during sampling. It can add or remove components through "birth" and "death" moves, using Bayesian inference to determine the most likely number of components.

### Mathematical Model

The target mixture distribution is defined as:

```
p(z) = Σₖ₌₁ᴷ pₖ · N(z | μₖ, σₖ²)
```

Where:

- N(z | μₖ, σₖ²) is the normal density for the k-th component
- pₖ represents the mixing proportion of the k-th component
- K is the number of components (inferred dynamically)
- μₖ are the component means
- σₖ are the component standard deviations

## 🛠️ Usage

1. **Upload Data**: Drag and drop a CSV file or click to browse
2. **Load Example**: Click "Load 3-Component Example" for a demonstration
3. **Configure Parameters**: Set MCMC iterations, burn-in, and thinning
4. **Run Analysis**: Click "Start Analysis" to begin RJMCMC sampling
5. **View Results**: Explore posterior distributions and convergence diagnostics

## 🔗 Related Links

- **[RJMC Package](https://github.com/dchodge/rjmc)** - Source R package
- **[Detailed Vignette](https://dchodge.github.io/rjmc/articles/Ex1_mixture.html)** - Technical documentation
- **[RJMC Documentation](https://dchodge.github.io/rjmc/)** - Package documentation

## 📄 License

This project is part of the RJMC package ecosystem. Please refer to the main package license.

## 🤝 Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

---

**Developed by David Hodgson**
*Part of the RJMC package for Reversible Jump MCMC analysis*
