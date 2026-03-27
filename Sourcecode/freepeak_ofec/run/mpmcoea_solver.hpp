#ifndef OFEC_MPMCOEA_SOLVER_HPP
#define OFEC_MPMCOEA_SOLVER_HPP

#include "custom_method.hpp" 
#include "../instance/problem/continuous/free_peaks/free_peaks.h"
#include "../core/environment/environment.h"
#include <queue>
#include <numeric>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>

namespace ofec {

    // NEAREST-BETTER CLUSTERING (NBC)

    std::vector<std::vector<int>> nearestBetterClustering(
        const std::vector<std::vector<double>>& population,
        const std::vector<double>& fitness,
        double phi = 2.0)
    {
        int n = static_cast<int>(population.size());
        if (n == 0) return {};

        // Build nearest-better graph
        struct Edge { int from, to; double weight; };
        std::vector<Edge> edges;

        for (int i = 0; i < n; i++) {
            double min_dist = std::numeric_limits<double>::max();
            int nb_idx = -1;
            for (int j = 0; j < n; j++) {
                if (fitness[j] > fitness[i] + 1e-9) {  // Strictly better (maximization)
                    double d = 0.0;
                    for (size_t k = 0; k < population[i].size(); k++) {
                        double diff = population[i][k] - population[j][k];
                        d += diff * diff;
                    }
                    d = std::sqrt(d);
                    if (d < min_dist) { min_dist = d; nb_idx = j; }
                }
            }
            if (nb_idx != -1) edges.push_back({ i, nb_idx, min_dist });
        }

        // Handle no better neighbors
        if (edges.empty()) {
            std::vector<std::vector<int>> clusters;
            for (int i = 0; i < n; i++) clusters.push_back({ i });
            return clusters;
        }

        // Prune edges
        double sum_weights = 0.0;
        for (const auto& e : edges) sum_weights += e.weight;
        double threshold = phi * (sum_weights / edges.size());

        // Build adjacency list
        std::vector<std::vector<int>> adj(n);
        for (const auto& e : edges) {
            if (e.weight <= threshold) {
                adj[e.from].push_back(e.to);
                adj[e.to].push_back(e.from);
            }
        }

        // BFS for connected components
        std::vector<std::vector<int>> clusters;
        std::vector<bool> visited(n, false);
        for (int i = 0; i < n; i++) {
            if (!visited[i]) {
                std::vector<int> cluster;
                std::queue<int> q;
                q.push(i); visited[i] = true;
                while (!q.empty()) {
                    int u = q.front(); q.pop();
                    cluster.push_back(u);
                    for (int v : adj[u]) {
                        if (!visited[v]) { visited[v] = true; q.push(v); }
                    }
                }
                clusters.push_back(cluster);
            }
        }
        return clusters;
    }

    // MPMMO BENCHMARK WRAPPER 

    class MPMMO_Benchmark {
    private:
        FreePeaks* m_problem;  
        Environment* m_env;
        int m_num_dms;
        std::vector<std::vector<int>> m_party_indices;

    public:
        MPMMO_Benchmark(FreePeaks* prob, int num_dms, Environment* env)
            : m_problem(prob), m_num_dms(num_dms), m_env(env) {
            decomposeVariables();
        }

        void decomposeVariables() {
            int total_vars = m_problem->numberVariables();
            m_party_indices.resize(m_num_dms);

            // For 2D problem with 2 DMs:
            // DM 0 gets variable 0 (x₁)
            // DM 1 gets variable 1 (x₂)

            for (int v = 0; v < total_vars; v++) {
                int dm_id = v % m_num_dms;  // Round-robin assignment
                m_party_indices[dm_id].push_back(v);
            }

            // DEBUG: Print variable assignment
            std::cout << "[DEBUG] Variable Decomposition:" << std::endl;
            for (int dm = 0; dm < m_num_dms; dm++) {
                std::cout << "  DM " << dm << " controls variables: [";
                for (size_t i = 0; i < m_party_indices[dm].size(); i++) {
                    std::cout << m_party_indices[dm][i];
                    if (i < m_party_indices[dm].size() - 1) std::cout << ", ";
                }
                std::cout << "]" << std::endl;
            }
        }

        // Context-based evaluation using PUBLIC FreePeaks API
        double evaluateDM(int dm_id, const std::vector<double>& my_vars,
            const std::vector<std::vector<double>>& context,
            Environment* env)
        {
            if (dm_id == 0 && my_vars.size() > 0) {
                std::cout << "[DEBUG evaluateDM] DM0 my_vars[0]=" << my_vars[0] << std::endl;
            }
            // Build full solution vector
            std::vector<double> full_sol(m_problem->numberVariables(), 0.5);

            const auto& my_indices = m_party_indices[dm_id];
            for (size_t i = 0; i < my_indices.size(); i++) {
                full_sol[my_indices[i]] = my_vars[i];
            }

            for (int other = 0; other < m_num_dms; other++) {
                if (other == dm_id) continue;
                const auto& other_indices = m_party_indices[other];
                for (size_t i = 0; i < other_indices.size(); i++) {
                    full_sol[other_indices[i]] = context[other][i];
                }
            }

            // Evaluate SAFELY via public API
            auto sol = m_problem->createSolution();
            auto& sol_vec = dynamic_cast<Solution<VariableVector<Real>>&>(*sol).variable().vector();
            sol_vec = full_sol;
            sol->evaluate(m_env, false);
            return sol->objective(0);   // Single objective (maximization)
        }

        int getNumDMs() const { return m_num_dms; }
        int getNumVariables() const { return m_problem->numberVariables(); }
        const std::vector<int>& getPartyIndices(int p) const { return m_party_indices[p]; }
        FreePeaks* getProblem() const { return m_problem; }
    };


    // MPM-CoEA SOLVER 

    class MPM_CoEA {
    private:
        std::shared_ptr<MPMMO_Benchmark> m_bench;
        std::shared_ptr<Environment> m_env;
        int m_pop_size, m_med_pop_size, m_K;

        // Competing populations: [DM][subpopulation_id][individual][variable]
        std::vector<std::vector<std::vector<std::vector<double>>>> m_competing_pops;
        std::vector<std::vector<std::vector<double>>> m_competing_fitness;

        // Mediating population: [individual][variable]
        std::vector<std::vector<double>> m_mediating_pop;
        std::vector<double> m_mediating_fitness;

        // Context vectors: Best representative per DM [DM][variable]
        std::vector<std::vector<double>> m_context_vectors;

        // Member variable
        std::string output_directory;  // For saving solutions

        


    public:

        // Constructor parameter:
        MPM_CoEA(std::shared_ptr<MPMMO_Benchmark> bench,
            std::shared_ptr<Environment> env,
            int pop_size = 50,
            int med_pop_size = 100,
            int K = 3,
            const std::string& output_dir = "visualization/solutions")
            : m_bench(bench), m_env(env), m_pop_size(pop_size),
            m_med_pop_size(med_pop_size), m_K(K), output_directory(output_dir)
        {
            // Creating output directory
            std::filesystem::create_directories(output_directory);
        }

        void saveSolutionsAtGeneration(int gen_num) {
            std::string filename = output_directory + "/gen_" +
                std::to_string(gen_num) + "_solutions.txt";

            std::ofstream outFile(filename);
            outFile << std::fixed << std::setprecision(6);

            // Write header
            outFile << "# Generation: " << gen_num << std::endl;
            outFile << "# Format: x1 x2 fitness dm_id" << std::endl;

            // Save mediating population solutions
            for (size_t i = 0; i < m_mediating_pop.size(); i++) {
                // Evaluate fitness
                auto sol = m_bench->getProblem()->createSolution();
                auto& sol_vec = dynamic_cast<Solution<VariableVector<Real>>&>(*sol).variable().vector();
                sol_vec[0] = static_cast<Real>(m_mediating_pop[i][0]);
                sol_vec[1] = static_cast<Real>(m_mediating_pop[i][1]);
                sol->evaluate(m_env.get(), false);

                outFile << m_mediating_pop[i][0] << " "
                    << m_mediating_pop[i][1] << " "
                    << sol->objective(0) << " "
                    << "mediating" << "\n";
            }

            outFile.close();
        }

        MPM_CoEA(std::shared_ptr<MPMMO_Benchmark> bench, std::shared_ptr<Environment> env,
            int pop_size = 50, int med_pop_size = 100, int K = 3)
            : m_bench(bench), m_env(env), m_pop_size(pop_size),
            m_med_pop_size(med_pop_size), m_K(K) {
        }

        void initialize() {
            int dim = m_bench->getNumVariables();
            const int NBC_POP_SIZE = 500;

            // Create large initial population for NBC
            std::vector<std::vector<double>> global_pop(NBC_POP_SIZE, std::vector<double>(dim));
            std::vector<std::vector<double>> global_fitness(m_bench->getNumDMs(),
                std::vector<double>(NBC_POP_SIZE));

            // Random initialization
            for (int i = 0; i < NBC_POP_SIZE; i++) {
                for (int d = 0; d < dim; d++) {
                    global_pop[i][d] = (double)rand() / RAND_MAX;
                }
                std::vector<std::vector<double>> neutral_context(m_bench->getNumDMs());
                for (int dm = 0; dm < m_bench->getNumDMs(); dm++) {
                    neutral_context[dm].resize(m_bench->getPartyIndices(dm).size(), 0.5);
                }
                
                for (int dm = 0; dm < m_bench->getNumDMs(); dm++) {
                    global_fitness[dm][i] = m_bench->evaluateDM(dm, global_pop[i], neutral_context, m_env.get());
                }
            }

            // Apply NBC per DM → create subpopulations
            m_competing_pops.resize(m_bench->getNumDMs());
            m_competing_fitness.resize(m_bench->getNumDMs());

            for (int dm = 0; dm < m_bench->getNumDMs(); dm++) {
                auto clusters = nearestBetterClustering(global_pop, global_fitness[dm], 2.0);
                m_competing_pops[dm].resize(clusters.size());
                m_competing_fitness[dm].resize(clusters.size());

                for (size_t c = 0; c < clusters.size(); c++) {
                    int take = std::min((int)clusters[c].size(), m_pop_size);
                    m_competing_pops[dm][c].resize(take, std::vector<double>(dim));
                    m_competing_fitness[dm][c].resize(take);

                    for (int i = 0; i < take; i++) {
                        m_competing_pops[dm][c][i] = global_pop[clusters[c][i]];
                        m_competing_fitness[dm][c][i] = global_fitness[dm][clusters[c][i]];
                    }
                }
            }

            // Initialize mediating population
            m_mediating_pop.resize(m_med_pop_size, std::vector<double>(dim));
            m_mediating_fitness.resize(m_med_pop_size);

            // Collect top candidates
            std::vector<std::pair<double, std::vector<double>>> candidates;
            for (int dm = 0; dm < m_bench->getNumDMs(); dm++) {
                for (size_t c = 0; c < m_competing_pops[dm].size(); c++) {
                    for (size_t i = 0; i < m_competing_pops[dm][c].size(); i++) {
                        candidates.emplace_back(m_competing_fitness[dm][c][i],
                            m_competing_pops[dm][c][i]);
                    }
                }
            }

            // Sort by fitness (descending for maximization)
            std::sort(candidates.begin(), candidates.end(),
                [](const auto& a, const auto& b) { return a.first > b.first; });

            // Select diverse top solutions
            int selected = 0;
            for (size_t i = 0; i < candidates.size() && selected < m_med_pop_size; i++) {
                bool is_duplicate = false;
                for (int j = 0; j < selected; j++) {
                    double dist = 0.0;
                    for (int d = 0; d < dim; d++) {
                        double diff = candidates[i].second[d] - m_mediating_pop[j][d];
                        dist += diff * diff;
                    }
                    if (std::sqrt(dist) < 0.05) { is_duplicate = true; break; }
                }
                if (!is_duplicate) {
                    m_mediating_pop[selected] = candidates[i].second;
                    m_mediating_fitness[selected] = candidates[i].first;
                    selected++;
                }
            }

            // Fill remaining slots randomly if needed
            while (selected < m_med_pop_size) {
                for (int d = 0; d < dim; d++) {
                    m_mediating_pop[selected][d] = (double)rand() / RAND_MAX;
                }
                std::vector<std::vector<double>> neutral_context(m_bench->getNumDMs(),
                    std::vector<double>(1, 0.5));
                m_mediating_fitness[selected] = m_bench->evaluateDM(0, m_mediating_pop[selected],
                    neutral_context, m_env.get());
                selected++;
            }

            // Initialize context vectors with RANDOM values (NOT 0.5!)
            m_context_vectors.resize(m_bench->getNumDMs());
            for (int dm = 0; dm < m_bench->getNumDMs(); dm++) {
                const auto& indices = m_bench->getPartyIndices(dm);
                m_context_vectors[dm].resize(indices.size());

                // Use RANDOM initialization, not fixed 0.5
                for (size_t i = 0; i < indices.size(); i++) {
                    m_context_vectors[dm][i] = (double)rand() / RAND_MAX;  // Random [0,1]
                }
            }
        }

        void run_one_generation() {
            static int gen_count = 0;

            // DEBUG: Print context at start of generation
            if (gen_count % 20 == 0) {
                std::cout << "\n[DEBUG gen " << gen_count << "] Context vectors:" << std::endl;
                for (int dm = 0; dm < m_bench->getNumDMs(); dm++) {
                    std::cout << "  DM" << dm << ": [";
                    for (size_t i = 0; i < m_context_vectors[dm].size(); i++) {
                        std::cout << std::fixed << std::setprecision(4)
                            << m_context_vectors[dm][i];
                        if (i < m_context_vectors[dm].size() - 1) std::cout << ", ";
                    }
                    std::cout << "]" << std::endl;
                }
            }
            gen_count++;

            int dim = m_bench->getNumVariables();
            int num_dms = m_bench->getNumDMs();

            // Step 1: Evolve competing populations per DM
            for (int dm = 0; dm < num_dms; dm++) {
                const auto& dm_indices = m_bench->getPartyIndices(dm);

                for (size_t subpop = 0; subpop < m_competing_pops[dm].size(); subpop++) {
                    // Simple Gaussian mutation on ALL variables
                    for (size_t i = 1; i < m_competing_pops[dm][subpop].size(); i++) {
                        for (int d = 0; d < dim; d++) {
                            double mutation_rate = (dm == 0) ? 0.3 : 0.5;  // Higher mutation for DM1
                            m_competing_pops[dm][subpop][i][d] += ((double)rand() / RAND_MAX - 0.5) * 0.3;

                            // Boundary handling [0,1]
                            if (m_competing_pops[dm][subpop][i][d] < 0)
                                m_competing_pops[dm][subpop][i][d] = 0;
                            if (m_competing_pops[dm][subpop][i][d] > 1)
                                m_competing_pops[dm][subpop][i][d] = 1;
                        }

                        // Evaluate with current context
                        m_competing_fitness[dm][subpop][i] =
                            m_bench->evaluateDM(dm, m_competing_pops[dm][subpop][i],
                                m_context_vectors, m_env.get());
                    }

                    // Sort subpopulation by fitness (elitism)
                    std::vector<size_t> indices(m_competing_pops[dm][subpop].size());
                    std::iota(indices.begin(), indices.end(), 0);
                    std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
                        return m_competing_fitness[dm][subpop][a] >
                            m_competing_fitness[dm][subpop][b];
                        });

                    // Reorder population
                    std::vector<std::vector<double>> new_pop;
                    std::vector<double> new_fit;
                    for (size_t idx : indices) {
                        new_pop.push_back(m_competing_pops[dm][subpop][idx]);
                        new_fit.push_back(m_competing_fitness[dm][subpop][idx]);
                    }
                    m_competing_pops[dm][subpop] = new_pop;
                    m_competing_fitness[dm][subpop] = new_fit;
                }
            }

            // Step 2: Update context vectors
            for (int dm = 0; dm < num_dms; dm++) {
                double best_fit = -1e9;
                std::vector<double> best_vars;
                const auto& indices = m_bench->getPartyIndices(dm);

                // Search ALL subpopulations for best individual
                for (size_t subpop = 0; subpop < m_competing_pops[dm].size(); subpop++) {
                    if (m_competing_pops[dm][subpop].size() > 0 &&
                        m_competing_fitness[dm][subpop][0] > best_fit)
                    {
                        best_fit = m_competing_fitness[dm][subpop][0];
                        best_vars.clear();

                        // Extract this DM's variables only from the full solution
                        for (size_t i = 0; i < indices.size(); i++) {
                            int var_idx = indices[i];
                            if (var_idx >= 0 && var_idx < dim) {
                                best_vars.push_back(m_competing_pops[dm][subpop][0][var_idx]);
                            }
                        }
                    }
                }

                // always update context
                if (best_vars.size() > 0) {
                    m_context_vectors[dm] = best_vars;
                }

                // DEBUG: Print context update
                if (gen_count % 20 == 0 && dm == 0) {
                    std::cout << "[DEBUG context update] DM" << dm << " context=[";
                    for (size_t i = 0; i < m_context_vectors[dm].size(); i++) {
                        std::cout << std::fixed << std::setprecision(4)
                            << m_context_vectors[dm][i];
                        if (i < m_context_vectors[dm].size() - 1) std::cout << ", ";
                    }
                    std::cout << "], fitness=" << best_fit << std::endl;
                }
            }

            // Step 3: Evolve mediating population with DE/rand/1/bin
            std::vector<std::vector<double>> new_med_pop = m_mediating_pop;
            std::vector<double> new_med_fit(m_med_pop_size);

            for (int i = 0; i < m_med_pop_size; i++) {
                // DE/rand/1/bin mutation
                int r1 = rand() % m_med_pop_size;
                int r2 = rand() % m_med_pop_size;
                int r3 = rand() % m_med_pop_size;
                while (r2 == r1) r2 = rand() % m_med_pop_size;
                while (r3 == r1 || r3 == r2) r3 = rand() % m_med_pop_size;

                double F = 0.5, CR = 0.9;
                for (int d = 0; d < dim; d++) {
                    if ((rand() / (double)RAND_MAX) < CR || d == dim - 1) {
                        new_med_pop[i][d] = m_mediating_pop[r1][d] +
                            F * (m_mediating_pop[r2][d] - m_mediating_pop[r3][d]);

                        // Boundary handling
                        if (new_med_pop[i][d] < 0) new_med_pop[i][d] = 0;
                        if (new_med_pop[i][d] > 1) new_med_pop[i][d] = 1;
                    }
                }

                // Use neutral_context so mediating pop's variables are actually evaluated
                //std::vector<std::vector<double>> neutral_context(num_dms, std::vector<double>(1, 0.5));
                new_med_fit[i] = m_bench->evaluateDM(0, new_med_pop[i],
                    m_context_vectors, m_env.get());  // Changed from m_context_vectors
            }

            // Selection
            for (int i = 0; i < m_med_pop_size; i++) {
                if (new_med_fit[i] > m_mediating_fitness[i]) {
                    m_mediating_pop[i] = new_med_pop[i];
                    m_mediating_fitness[i] = new_med_fit[i];
                }
            }
        }

        double getBestConsensusFitness() const {
            if (m_mediating_fitness.empty()) return -1e9;
            return *std::max_element(m_mediating_fitness.begin(),
                m_mediating_fitness.end());
        }

        const std::vector<std::vector<double>>& getMediatingPopulation() const {
            return m_mediating_pop;
        }
    };


    // CEC2015 METRICS EVALUATION 

    struct CEC2015Metrics {
        double SR = 0.0;    // Success Rate (%)
        double ANOF = 0.0;  // Average Number of Optima Found
        double SP = 0.0;    // Success Performance
        double MPR = 0.0;   // Maximum Peak Ratio
    };

    

    // CEC2015 METRICS EVALUATION (Complete Corrected Version)
    CEC2015Metrics evaluateCEC2015Metrics(
        const std::vector<std::vector<double>>& candidates,
        FreePeaks* problem,
        ofec::Environment* env,
        double epsilon_x = 1e-2,   
        double epsilon_f = 1.0)   
    {
        CEC2015Metrics metrics;

        // Get optima from problem
        auto optima_ptr = problem->optima();

        // DEBUG: Print optima information
        std::cout << "\n>>> [DEBUG] Optima Information:" << std::endl;
        if (!optima_ptr) {
            std::cout << "    [WARNING] problem->optima() returned nullptr!" << std::endl;
            std::cout << "    This means optima were not calculated during bindData()" << std::endl;
            return metrics; // Returns zeros
        }

        int num_optima = optima_ptr->numberSolutions();
        std::cout << "    Total optima found: " << num_optima << std::endl;

        if (num_optima == 0) {
            std::cout << "    [ERROR] No optima stored in problem!" << std::endl;
            std::cout << "    Check if problem->bindData() was called after loading mpmmo_1.txt" << std::endl;
            return metrics;
        }

        // Print optima positions and fitness
        for (int i = 0; i < num_optima && i < 5; i++) { // Print first 5 optima
            auto& opt_sol = optima_ptr->solution(i);
            auto& opt_vars = opt_sol.variable().vector();
            double opt_fitness = opt_sol.objective(0);

            std::cout << "    Optimum " << i << ": fitness=" << std::fixed << std::setprecision(4)
                << opt_fitness << ", position=[";
            for (size_t d = 0; d < opt_vars.size(); d++) {
                std::cout << opt_vars[d];
                if (d < opt_vars.size() - 1) std::cout << ", ";
            }
            std::cout << "]" << std::endl;
        }
        if (num_optima > 5) {
            std::cout << "    ... and " << (num_optima - 5) << " more optima" << std::endl;
        }

        // Track found optima
        std::vector<bool> optima_found(num_optima, false);
        int total_found = 0;

        // DEBUG: Print candidate information
        std::cout << "\n>>> [DEBUG] Candidate Solutions:" << std::endl;
        std::cout << "    Total candidates: " << candidates.size() << std::endl;

        for (size_t c = 0; c < candidates.size() && c < 5; c++) { // Print first 5 candidates
            // Evaluate candidate fitness
            auto sol = problem->createSolution();
            auto& sol_typed = dynamic_cast<ofec::Solution<ofec::VariableVector<ofec::Real>>&>(*sol);
            sol_typed.variable().vector() = candidates[c];
            sol->evaluate(env, false);
            double cand_fitness = sol->objective(0);

            std::cout << "    Candidate " << c << ": fitness=" << std::fixed << std::setprecision(4)
                << cand_fitness << ", position=[";
            for (size_t d = 0; d < candidates[c].size(); d++) {
                std::cout << candidates[c][d];
                if (d < candidates[c].size() - 1) std::cout << ", ";
            }
            std::cout << "]" << std::endl;

            // Check against all true optima
            for (int i = 0; i < num_optima; i++) {
                if (optima_found[i]) continue;

                // Position distance
                auto& opt_vars = optima_ptr->solution(i).variable().vector();
                double pos_dist = 0.0;
                for (size_t d = 0; d < candidates[c].size(); d++) {
                    double diff = candidates[c][d] - opt_vars[d];
                    pos_dist += diff * diff;
                }
                pos_dist = std::sqrt(pos_dist);

                // Fitness distance
                double opt_fitness = optima_ptr->solution(i).objective(0);
                double fit_dist = std::abs(cand_fitness - opt_fitness);

                // Check if within tolerance
                if (pos_dist < epsilon_x && fit_dist < epsilon_f) {
                    optima_found[i] = true;
                    total_found++;
                    std::cout << "        -> MATCHED optimum " << i << " (pos_dist=" << std::setprecision(6)
                        << pos_dist << ", fit_dist=" << fit_dist << ")" << std::endl;
                    break;
                }
            }
        }
        if (candidates.size() > 5) {
            std::cout << "    ... and " << (candidates.size() - 5) << " more candidates" << std::endl;
        }

        std::cout << "\n>>> [DEBUG] Matching Summary:" << std::endl;
        std::cout << "    Optima found: " << total_found << "/" << num_optima << std::endl;
        std::cout << "    Tolerance: epsilon_x=" << epsilon_x << ", epsilon_f=" << epsilon_f << std::endl;

        // Compute metrics
        metrics.ANOF = static_cast<double>(total_found) / num_optima;
        metrics.SR = (total_found == num_optima) ? 100.0 : 0.0;
        metrics.SP = (metrics.SR > 0) ? (2000.0 * problem->numberVariables() * num_optima / 100.0)
            : std::numeric_limits<double>::infinity();

        // MPR calculation (MAXIMIZATION version)
        double sum_num = 0.0, sum_den = 0.0;
        double global_opt_fitness = -std::numeric_limits<double>::max();

        for (int i = 0; i < num_optima; i++) {
            double f = optima_ptr->solution(i).objective(0);
            if (f > global_opt_fitness) global_opt_fitness = f;
        }

        for (int i = 0; i < num_optima; i++) {
            double true_fitness = optima_ptr->solution(i).objective(0);
            double val = (global_opt_fitness - true_fitness) + 1.0;
            if (optima_found[i]) {
                sum_num += val;
            }
            sum_den += val;
        }
        metrics.MPR = (sum_den > 0) ? (sum_num / sum_den) : 0.0;

        return metrics;
    }

    // CREATE SIMPLE 2D MPMMO PROBLEM FOR VISUALIZATION

    void createSimple2DProblem(const std::string& problem_name = "mpmmo_simple_2d") {
        using namespace ofec;
        using namespace free_peaks;

        std::string dirname = "multiparty_multimodal/";
        int numDim = 2, numObj = 1, numCon = 0;

        std::shared_ptr<ofec::Random> rnd(new Random(0.42));

        FreePeaks::registerFP();
        std::string freepeakName = "free_peaks";
        std::shared_ptr<Environment> env(generateEnvironmentByFactory(freepeakName));
        env->recordInputParameters();
        env->initialize();
        env->setProblem(generateProblemByFactory(freepeakName));

        auto freepeak = CAST_FPs(env->problem());

        ParameterMap freepeak_param;
        freepeak_param["generation_type"] = std::string("assigned");
        freepeak_param["dataFile1"] = dirname + "/" + problem_name + ".txt";
        freepeak->inputParameters().input(freepeak_param);
        freepeak->initialize(env.get());

        std::shared_ptr<ofec::Random> pro_rnd(new Random(0.42));
        freepeak->setRandom(pro_rnd);
        freepeak->setSizes(numDim, numObj, numCon);

        // KD-TREE: Defines the search space decomposition
        // This defines the multimodal landscape (multiple optima)
        freepeak->setKDtree({
            {"root", {
                {"subspace_1", 0.5},
                {"subspace_2", 0.5}
            }}
            });

        // Create subproblems for the KD-tree subspaces
        // BOTH subspaces represent the SAME objective function (multi-modal)
        std::vector<std::string> subspaces = { "subspace_1", "subspace_2" };

        // Define MULTIPLE OPTIMA for the multimodal landscape
        // These are the optima that BOTH parties need to find together
        std::vector<std::vector<double>> peak_positions = {
            {0.3, 0.7},  // Global optimum at [0.3, 0.7]
            {0.7, 0.3},  // Local optimum at [0.7, 0.3]
            {0.2, 0.2},  // Additional local optimum
            {0.8, 0.8}   // Additional local optimum
        };
        std::vector<double> peak_heights = { 90.0, 70.0, 50.0, 60.0 };

        // Create a subproblem for each subspace in the KD-tree
        for (size_t s = 0; s < subspaces.size(); s++) {
            std::string subspace_name = subspaces[s];
            ParameterMap subpro_param;
            subpro_param["subspace"] = subspace_name;
            subpro_param["generation_type"] = std::string("assigned");
            subpro_param["dataFile1"] = dirname + "/" + subspace_name + ".txt";

            auto subpro(Subproblem::create());
            subpro->initialize(subpro_param, freepeak);

            // Distance: Euclidean
            {
                auto dis(FactoryFP<DistanceBase>::produce("Euclidean"));
                ParameterMap dis_param;
                dis->initialize(freepeak, subspace_name, dis_param);
                subpro->setDistance(dis);
            }

            // Function: One Peak Function
            // Each subspace gets ALL peaks (multi-modal)
            {
                ParameterMap fun_param;
                fun_param["generation_type"] = std::string("assigned");
                fun_param["dataFile1"] = dirname + "/" + subspace_name + "_onepeak.txt";

                auto func(FactoryFP<FunctionBase>::produce("one_peak"));
                func->initialize(freepeak, subspace_name, fun_param);
                auto onepeak_func = dynamic_cast<ofec::free_peaks::OnePeakFunction*>(func);

                // Add ALL peaks to this subspace (multi-modal)
                for (size_t p = 0; p < peak_positions.size(); p++) {
                    auto onepeak(FactoryFP<OnePeakBase>::produce("s1"));
                    ParameterMap onepeak_param;
                    onepeak_param["center_type"] = std::string("assigned");
                    onepeak_param["height"] = peak_heights[p];

                    std::vector<Real> peak_position = {
                        static_cast<Real>(peak_positions[p][0]),
                        static_cast<Real>(peak_positions[p][1])
                    };
                    onepeak_param["center_postion"] = peak_position;

                    onepeak->initialize(freepeak, subspace_name, onepeak_param);
                    onepeak_func->addOnePeaks(onepeak);
                }

                subpro->setFunction(func);
            }

            freepeak->setSubproblem(subspace_name, subpro);
        }

        freepeak->bindData();

        // Create directories and save files
        std::cout << "freepeak path\t" << FreePeaks::directory() + dirname << std::endl;
        std::filesystem::create_directories(FreePeaks::directory() + dirname);
        std::cout << "Subproblem path\t" << Subproblem::directory() + dirname << std::endl;
        std::filesystem::create_directories(Subproblem::directory() + dirname);
        std::cout << "OnePeakFunction path\t" << OnePeakFunction::directory() + dirname << std::endl;
        std::filesystem::create_directories(OnePeakFunction::directory() + dirname);

        // Set generation type to read_file to save parameters
        freepeak->inputParameters().at("generation_type")->setValue("read_file");
        for (auto& it : freepeak->subspaceTree().name_box_subproblem) {
            if (it.second.second != nullptr) {
                it.second.second->inputParameters().at("generation_type")->setValue("read_file");
                it.second.second->function()->inputParameters().at("generation_type")->setValue("read_file");
            }
        }

        freepeak->recordInputParameters();
        freepeak->outputTotalFile();

        std::cout << ">>> Simple 2D MPMMO problem created: " << problem_name << std::endl;
        std::cout << "    Number of optima: " << peak_positions.size() << std::endl;
        for (size_t i = 0; i < peak_positions.size(); i++) {
            std::cout << "    Optimum " << i << ": ["
                << peak_positions[i][0] << ", "
                << peak_positions[i][1] << "] fitness="
                << peak_heights[i] << std::endl;
        }
    }

    // OUTPUT FITNESS LANDSCAPE FOR VISUALIZATION
    void outputFitnessLandscape(ofec::Environment* env, const std::string& filepath, int resolution = 100) {
        using namespace ofec;

        auto con_pro = CAST_CONOP(env->problem());
        std::vector<std::vector<double>> fitness_grid(resolution, std::vector<double>(resolution));

        std::cout << "\n>>> Generating fitness landscape (" << resolution << "x" << resolution << " grid)..." << std::endl;

        for (int i = 0; i < resolution; i++) {
            for (int j = 0; j < resolution; j++) {
                // Create solution at grid point
                auto sol = env->problem()->createSolution();
                auto& sol_vec = dynamic_cast<Solution<VariableVector<Real>>&>(*sol).variable().vector();

                // Map grid coordinates to variable range [0, 1]
                sol_vec[0] = static_cast<Real>(i) / (resolution - 1);
                sol_vec[1] = static_cast<Real>(j) / (resolution - 1);

                // Evaluate
                sol->evaluate(env, false);
                fitness_grid[i][j] = sol->objective(0);
            }

            if (i % 10 == 0) {
                std::cout << "    Progress: " << (i * 100 / resolution) << "%" << std::endl;
            }
        }

        // Save to file
        std::ofstream outFile(filepath);
        outFile << std::fixed << std::setprecision(6);

        // Write data
        for (int i = 0; i < resolution; i++) {
            for (int j = 0; j < resolution; j++) {
                double x1 = static_cast<double>(i) / (resolution - 1);
                double x2 = static_cast<double>(j) / (resolution - 1);
                outFile << x1 << " " << x2 << " " << fitness_grid[i][j] << "\n";
            }
        }

        outFile.close();
        std::cout << ">>> Landscape saved to: " << filepath << std::endl;
    }

    // MAIN SOLVER ENTRY POINT WITH VISUALIZATION
    void runMPMCoEAExperiment() {
        srand(static_cast<unsigned int>(time(nullptr)));

        ofec::registerInstance();

        std::cout << "\n>>> [MPM-CoEA] Starting experiment with visualization" << std::endl;

        // Set working directory
        ofec::g_working_directory = R"(E:\HITSZ\Research\Multimodal_Multiparty_Optimization\ThesisProject\Data\ofec_data_new)";

        // Define output directory for visualization
        std::string output_dir = "E:/HITSZ/Research/Multimodal_Multiparty_Optimization/ThesisProject/Visualization";
        std::filesystem::create_directories(output_dir);

        // Step 0: Create simple 2D problem for visualization
        std::cout << "\n>>> Creating simple 2D problem..." << std::endl;
        createSimple2DProblem("simple_2d");

        // Step 1: Load the problem
        std::string freepeakName = "free_peaks";

        auto env = generateEnvironmentByFactory(freepeakName);
        env->recordInputParameters();
        env->initialize();
        env->setProblem(generateProblemByFactory(freepeakName));

        ParameterMap params;
        params["generation_type"] = std::string("read_file");
        params["dataFile1"] = std::string("visualization/simple_2d.txt");
        env->problem()->inputParameters().input(params);
        env->problem()->recordInputParameters();
        env->initializeProblem(0.5);

        auto free_peaks_ptr = CAST_FPs(env->problem());
        if (!free_peaks_ptr) {
            throw std::runtime_error("Failed to cast problem to FreePeaks");
        }

        std::cout << ">>> [MPM-CoEA] Loaded benchmark: simple_2d" << std::endl;
        std::cout << "    Dimensions: " << free_peaks_ptr->numberVariables()
            << ", Decision-makers: 2" << std::endl;

        // Step 2: Output fitness landscape (BEFORE optimization)
        std::cout << "\n>>> Generating fitness landscape..." << std::endl;
        outputFitnessLandscape(env, output_dir + "/landscape_before.txt", 100);  

        // Step 3: Create benchmark wrapper and solver
        auto benchmark = std::make_shared<MPMMO_Benchmark>(free_peaks_ptr, 2, env);
        auto env_shared = std::shared_ptr<Environment>(env, [](Environment*) {});  // No-op deleter
        MPM_CoEA solver(benchmark, env_shared, 100, 300, 5, output_dir + "/solutions"); 

        solver.initialize();

        // Step 4: Save initial population (generation 0)
        solver.saveSolutionsAtGeneration(0);

        // Step 5: Run optimization with solution tracking
        int max_generations = 200;
        std::cout << "\n>>> [MPM-CoEA] Starting optimization (" << max_generations << " generations)" << std::endl;
        for (int gen = 0; gen < max_generations; gen++) {
            solver.run_one_generation();

            if (gen % 10 == 0 || gen == max_generations - 1) {
                std::cout << "Gen " << gen + 1 << " | Best Consensus Fitness: "
                    << std::fixed << std::setprecision(4)
                    << solver.getBestConsensusFitness() << std::endl;
            }

            if (gen % 5 == 0 || gen == max_generations - 1) {
                solver.saveSolutionsAtGeneration(gen + 1);
            }
        }

        // Step 6: Output final landscape (AFTER optimization)
        std::cout << "\n>>> Generating final fitness landscape..." << std::endl;
        outputFitnessLandscape(env, output_dir + "/landscape_after.txt", 100);  

        std::cout << "\n>>> [MPM-CoEA] Experiment complete!" << std::endl;
        std::cout << "    Visualization files saved to:" << std::endl;
        std::cout << "    - " << output_dir << std::endl;
    }

    

} // namespace ofec

#endif // !OFEC_MPMCOEA_SOLVER_HPP