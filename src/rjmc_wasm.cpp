#include <emscripten.h>
#include <emscripten/bind.h>
#include <vector>
#include <string>
#include <iostream>
#include <random>
#include <cmath>
#include <algorithm>

// Simplified RJMCMC implementation for WebAssembly
class RJMCMC_WASM {
private:
    std::vector<double> data;
    std::vector<std::vector<double>> jump_matrix; // [p, mu, sigma] for each component
    std::vector<double> fitted_params;
    int n_components;
    int n_iterations;
    int burnin;
    int thin;
    double current_log_posterior;
    std::mt19937 rng;
    
    // Model parameters
    double mu_min = 0.0, mu_max = 20.0;
    double sigma_min = 0.3, sigma_max = 3.0;
    int max_components = 8;
    
public:
    RJMCMC_WASM() : rng(std::random_device{}()) {}
    
    void setData(const std::vector<double>& input_data) {
        data = input_data;
    }
    
    void setParameters(int iterations, int burnin_period, int thinning) {
        n_iterations = iterations;
        burnin = burnin_period;
        thin = thinning;
    }
    
    double uniform(double min_val, double max_val) {
        std::uniform_real_distribution<double> dist(min_val, max_val);
        return dist(rng);
    }
    
    int uniform_int(int min_val, int max_val) {
        std::uniform_int_distribution<int> dist(min_val, max_val);
        return dist(rng);
    }
    
    double normal(double mean, double std) {
        std::normal_distribution<double> dist(mean, std);
        return dist(rng);
    }
    
    double log_likelihood() {
        double ll = 0.0;
        for (double obs : data) {
            double mixture_prob = 0.0;
            for (int k = 0; k < n_components; k++) {
                double p = jump_matrix[0][k];
                double mu = jump_matrix[1][k];
                double sigma = jump_matrix[2][k];
                mixture_prob += p * (1.0 / (sigma * sqrt(2.0 * M_PI))) * 
                               exp(-0.5 * pow((obs - mu) / sigma, 2));
            }
            ll += log(mixture_prob);
        }
        return ll;
    }
    
    double log_prior() {
        double lp = 0.0;
        
        // Prior for number of components
        lp += log(1.0 / max_components);
        
        // Priors for mixture parameters
        for (int k = 0; k < n_components; k++) {
            double p = jump_matrix[0][k];
            double mu = jump_matrix[1][k];
            double sigma = jump_matrix[2][k];
            
            // Uniform priors
            if (p >= 0.0 && p <= 1.0) lp += 0.0; // log(1) = 0
            else lp = -std::numeric_limits<double>::infinity();
            
            if (mu >= mu_min && mu <= mu_max) lp += 0.0;
            else lp = -std::numeric_limits<double>::infinity();
            
            if (sigma >= sigma_min && sigma <= sigma_max) lp += 0.0;
            else lp = -std::numeric_limits<double>::infinity();
        }
        
        return lp;
    }
    
    void initialize() {
        // Initialize with 2 components
        n_components = 2;
        jump_matrix.resize(3);
        jump_matrix[0] = {0.5, 0.5}; // mixing proportions
        jump_matrix[1] = {1.0, 4.0}; // means
        jump_matrix[2] = {1.0, 1.0}; // standard deviations
        
        fitted_params = {0.0}; // sigma parameter
        current_log_posterior = log_likelihood() + log_prior();
    }
    
    void birth_proposal() {
        if (n_components >= max_components) return;
        
        // Sample new component parameters
        double p_new = uniform(0.0, 1.0);
        double mu_new = uniform(mu_min, mu_max);
        double sigma_new = uniform(sigma_min, sigma_max);
        
        // Adjust existing proportions
        for (int k = 0; k < n_components; k++) {
            jump_matrix[0][k] *= (1.0 - p_new);
        }
        
        // Add new component
        jump_matrix[0].push_back(p_new);
        jump_matrix[1].push_back(mu_new);
        jump_matrix[2].push_back(sigma_new);
        
        n_components++;
    }
    
    void death_proposal() {
        if (n_components <= 1) return;
        
        int idx = uniform_int(0, n_components - 1);
        double p_removed = jump_matrix[0][idx];
        
        // Remove component
        jump_matrix[0].erase(jump_matrix[0].begin() + idx);
        jump_matrix[1].erase(jump_matrix[1].begin() + idx);
        jump_matrix[2].erase(jump_matrix[2].begin() + idx);
        
        // Renormalize proportions
        for (int k = 0; k < n_components - 1; k++) {
            jump_matrix[0][k] /= (1.0 - p_removed);
        }
        
        n_components--;
    }
    
    void jump_proposal() {
        int idx = uniform_int(0, n_components - 1);
        
        // Update mixing proportion
        double p_new = jump_matrix[0][idx] + normal(0.0, 0.01);
        p_new = std::max(0.0, std::min(1.0, p_new));
        
        double diff = (jump_matrix[0][idx] - p_new) / (n_components - 1);
        for (int k = 0; k < n_components; k++) {
            if (k == idx) {
                jump_matrix[0][k] = p_new;
            } else {
                jump_matrix[0][k] += diff;
            }
        }
        
        // Update mean and sigma
        jump_matrix[1][idx] += normal(0.0, 0.1);
        jump_matrix[2][idx] = std::max(sigma_min, jump_matrix[2][idx] + normal(0.0, 0.1));
    }
    
    std::vector<std::vector<double>> run_mcmc() {
        initialize();
        
        std::vector<std::vector<double>> results;
        std::vector<double> component_counts;
        
        for (int iter = 0; iter < n_iterations; iter++) {
            // Propose move
            double move_type = uniform(0.0, 1.0);
            std::vector<std::vector<double>> old_jump = jump_matrix;
            int old_n_components = n_components;
            
            if (n_components == 1) {
                if (move_type < 0.5) {
                    birth_proposal();
                } else {
                    jump_proposal();
                }
            } else if (n_components == max_components) {
                if (move_type < 0.5) {
                    death_proposal();
                } else {
                    jump_proposal();
                }
            } else {
                if (move_type < 0.33) {
                    death_proposal();
                } else if (move_type < 0.67) {
                    jump_proposal();
                } else {
                    birth_proposal();
                }
            }
            
            // Calculate acceptance probability
            double new_log_posterior = log_likelihood() + log_prior();
            double alpha = std::min(1.0, exp(new_log_posterior - current_log_posterior));
            
            // Accept or reject
            if (uniform(0.0, 1.0) < alpha) {
                current_log_posterior = new_log_posterior;
            } else {
                jump_matrix = old_jump;
                n_components = old_n_components;
            }
            
            // Store results after burnin
            if (iter >= burnin && iter % thin == 0) {
                std::vector<double> sample;
                sample.push_back(n_components);
                for (int k = 0; k < n_components; k++) {
                    sample.push_back(jump_matrix[0][k]); // p
                    sample.push_back(jump_matrix[1][k]); // mu
                    sample.push_back(jump_matrix[2][k]); // sigma
                }
                results.push_back(sample);
                component_counts.push_back(n_components);
            }
        }
        
        return results;
    }
    
    std::vector<double> get_component_counts() {
        std::vector<double> counts(max_components + 1, 0.0);
        auto results = run_mcmc();
        
        for (const auto& sample : results) {
            int n_comp = static_cast<int>(sample[0]);
            if (n_comp <= max_components) {
                counts[n_comp]++;
            }
        }
        
        // Normalize
        double total = results.size();
        for (double& count : counts) {
            count /= total;
        }
        
        return counts;
    }
};

// Global instance
RJMCMC_WASM rjmc_instance;

extern "C" {
    EMSCRIPTEN_KEEPALIVE
    void set_data(double* data_ptr, int n_data) {
        std::vector<double> data_vec(data_ptr, data_ptr + n_data);
        rjmc_instance.setData(data_vec);
    }
    
    EMSCRIPTEN_KEEPALIVE
    void set_parameters(int iterations, int burnin, int thin) {
        rjmc_instance.setParameters(iterations, burnin, thin);
    }
    
    EMSCRIPTEN_KEEPALIVE
    double* run_mcmc_analysis() {
        auto results = rjmc_instance.run_mcmc();
        
        // Flatten results for JavaScript
        static std::vector<double> flat_results;
        flat_results.clear();
        
        for (const auto& sample : results) {
            for (double val : sample) {
                flat_results.push_back(val);
            }
        }
        
        return flat_results.data();
    }
    
    EMSCRIPTEN_KEEPALIVE
    int get_results_size() {
        auto results = rjmc_instance.run_mcmc();
        int total_size = 0;
        for (const auto& sample : results) {
            total_size += sample.size();
        }
        return total_size;
    }
    
    EMSCRIPTEN_KEEPALIVE
    double* get_component_counts() {
        static std::vector<double> counts = rjmc_instance.get_component_counts();
        return counts.data();
    }
}

// Emscripten bindings for easier JavaScript integration
using namespace emscripten;

EMSCRIPTEN_BINDINGS(rjmc_module) {
    class_<RJMCMC_WASM>("RJMCMC")
        .constructor<>()
        .function("setData", &RJMCMC_WASM::setData)
        .function("setParameters", &RJMCMC_WASM::setParameters)
        .function("runMCMC", &RJMCMC_WASM::run_mcmc)
        .function("getComponentCounts", &RJMCMC_WASM::get_component_counts);
}
