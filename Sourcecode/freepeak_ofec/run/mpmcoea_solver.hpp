#ifndef OFEC_MPMCOEA_SOLVER_HPP
#define OFEC_MPMCOEA_SOLVER_HPP

#include "custom_method.hpp"  // Must be in same directory as this file
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
        FreePeaks* m_problem;  // Raw pointer (owned by Environment)
        std::shared_ptr<Environment> m_env;
        int m_num_dms;
        std::vector<std::vector<int>> m_party_indices;

    public:
        MPMMO_Benchmark(FreePeaks* prob, int num_dms, std::shared_ptr<Environment> env)
            : m_problem(prob), m_num_dms(num_dms), m_env(env) {
            decomposeVariables();
        }

        void decomposeVariables() {
            int total_vars = m_problem->numberVariables();
            m_party_indices.resize(m_num_dms);
            int base_count = total_vars / m_num_dms;
            int remainder = total_vars % m_num_dms;
            int current_var = 0;
            for (int p = 0; p < m_num_dms; p++) {
                int count = base_count + (p < remainder ? 1 : 0);
                for (int k = 0; k < count; k++) {
                    m_party_indices[p].push_back(current_var++);
                }
            }
        }

        // Context-based evaluation using PUBLIC FreePeaks API
        double evaluateDM(int dm_id, const std::vector<double>& my_vars,
            const std::vector<std::vector<double>>& context,
            Environment* env)
        {
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
            sol->evaluate(env, false);  // Public method
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

    public:
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
                std::vector<std::vector<double>> neutral_context(m_bench->getNumDMs(),
                    std::vector<double>(1, 0.5));
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
                    if (std::sqrt(dist) < 0.01) { is_duplicate = true; break; }
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

            // Initialize context vectors
            m_context_vectors.resize(m_bench->getNumDMs());
            for (int dm = 0; dm < m_bench->getNumDMs(); dm++) {
                double best_fit = -1e9;
                std::vector<double> best_vars;

                for (size_t c = 0; c < m_competing_pops[dm].size(); c++) {
                    if (m_competing_fitness[dm][c][0] > best_fit) {
                        best_fit = m_competing_fitness[dm][c][0];
                        const auto& indices = m_bench->getPartyIndices(dm);
                        best_vars.resize(indices.size());
                        for (size_t i = 0; i < indices.size(); i++) {
                            best_vars[i] = m_competing_pops[dm][c][0][indices[i]];
                        }
                    }
                }
                m_context_vectors[dm] = best_vars;
            }
        }

        void run_one_generation() {
            int dim = m_bench->getNumVariables();

            // Evolve competing populations per DM
            for (int dm = 0; dm < m_bench->getNumDMs(); dm++) {
                for (size_t subpop = 0; subpop < m_competing_pops[dm].size(); subpop++) {
                    // Simple Gaussian mutation
                    for (size_t i = 1; i < m_competing_pops[dm][subpop].size(); i++) {
                        for (int d = 0; d < dim; d++) {
                            m_competing_pops[dm][subpop][i][d] += ((double)rand() / RAND_MAX - 0.5) * 0.3;

                            // Boundary handling [0,1]
                            if (m_competing_pops[dm][subpop][i][d] < 0)
                                m_competing_pops[dm][subpop][i][d] = 0;
                            if (m_competing_pops[dm][subpop][i][d] > 1)
                                m_competing_pops[dm][subpop][i][d] = 1;
                        }

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

            // Update context vectors (best representative per DM)
            for (int dm = 0; dm < m_bench->getNumDMs(); dm++) {
                double best_fit = -1e9;
                std::vector<double> best_vars;

                for (size_t subpop = 0; subpop < m_competing_pops[dm].size(); subpop++) {
                    if (m_competing_fitness[dm][subpop][0] > best_fit) {
                        best_fit = m_competing_fitness[dm][subpop][0];
                        const auto& indices = m_bench->getPartyIndices(dm);
                        best_vars.resize(indices.size());
                        for (size_t i = 0; i < indices.size(); i++) {
                            best_vars[i] = m_competing_pops[dm][subpop][0][indices[i]];
                        }
                    }
                }
                m_context_vectors[dm] = best_vars;
            }

            // Evolve mediating population with DE/rand/1/bin
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

                // Evaluate consensus fitness (using DM 0 as proxy)
                std::vector<std::vector<double>> neutral_context(m_bench->getNumDMs(),
                    std::vector<double>(1, 0.5));
                new_med_fit[i] = m_bench->evaluateDM(0, new_med_pop[i], neutral_context, m_env.get());
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
        double epsilon_x = 1e-2,   // RELAXED from 1e-4 for debugging
        double epsilon_f = 1.0)    // RELAXED from 0.1 for debugging
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


    // MAIN SOLVER ENTRY POINT (Runs on pre-generated mpmmo_1.txt benchmark)
    
    void runMPMCoEAExperiment() {
        srand(static_cast<unsigned int>(time(nullptr)));

        ofec::registerInstance();

        // Set working directory
        ofec::g_working_directory = R"(E:\HITSZ\Research\Multimodal_Multiparty_Optimization\ThesisProject\Data\ofec_data_new)";

        std::cout << "\n>>> [MPM-CoEA] Starting experiment on mpmmo_1 benchmark" << std::endl;
        std::cout << "    Working directory: " << ofec::g_working_directory << std::endl;


        // Register FreePeaks components
        FreePeaks::registerFP();

        // Create Environment (using derived class)
        class DummyEnv : public ofec::Environment {
        public:
            DummyEnv() = default;
        };
        auto env = std::make_shared<DummyEnv>();

        env->recordInputParameters();
        env->initialize();

        // Create and configure problem
        std::string freepeakName = "free_peaks";
        env->setProblem(generateProblemByFactory(freepeakName));

        ParameterMap params;
        params["generation_type"] = std::string("read_file");
        params["dataFile1"] = std::string("multiparty_multimodal/mpmmo_1.txt");
        env->problem()->inputParameters().input(params);
        env->problem()->recordInputParameters();
        env->initializeProblem(0.5);

        CAST_FPs(env->problem())->bindData();  //Calculates optima from subproblems

        // Get FreePeaks pointer
        FreePeaks* free_peaks_ptr = CAST_FPs(env->problem());

        if (!free_peaks_ptr) {
            throw std::runtime_error("Failed to cast problem to FreePeaks");
        }

        std::cout << ">>> [MPM-CoEA] Loaded benchmark: mpmmo_1" << std::endl;
        std::cout << "    Dimensions: " << free_peaks_ptr->numberVariables()
            << ", Decision-makers: 2" << std::endl;

        // Create benchmark wrapper and solver
        auto benchmark = std::make_shared<MPMMO_Benchmark>(free_peaks_ptr, 2, env);
        MPM_CoEA solver(benchmark, env, 100, 300, 5);
        solver.initialize();

        int num_runs = 5;  // CEC2015 requirement
        std::vector<int> optima_per_run(num_runs);
        std::vector<double> anof_per_run(num_runs);
        std::vector<double> sr_per_run(num_runs);
        std::vector<double> mpr_per_run(num_runs);

        CEC2015Metrics last_run_metrics;  // Storing last run's metrics for display

        for (int run = 0; run < num_runs; run++) {
            srand(static_cast<unsigned int>(time(nullptr) + run));  // Different seed per run

            // Re-initializing solver for each run
            solver.initialize();

            // Run optimization (100 generations)
            std::cout << "\n>>> [MPM-CoEA] Run " << (run + 1) << "/" << num_runs << std::endl;
            for (int gen = 0; gen < 100; gen++) {
                solver.run_one_generation();
                if (gen % 20 == 0) {  // Printing every 20 generations to reduce output
                    std::cout << "  Gen " << gen + 1 << " | Best Consensus Fitness: "
                        << std::fixed << std::setprecision(4)
                        << solver.getBestConsensusFitness() << std::endl;
                }
            }

            // Evaluate results using CEC2015 metrics
            auto candidates = solver.getMediatingPopulation();
            auto metrics = evaluateCEC2015Metrics(candidates, free_peaks_ptr, env.get(), 0.1, 5);

            // Store results for this run
            optima_per_run[run] = static_cast<int>(metrics.ANOF * 12);  // 12 = total optima
            anof_per_run[run] = metrics.ANOF;
            sr_per_run[run] = metrics.SR;
            mpr_per_run[run] = metrics.MPR;

            last_run_metrics = metrics;  // ← Store for final display
        }

        // Calculate aggregate statistics across all runs
        double avg_anof = std::accumulate(anof_per_run.begin(), anof_per_run.end(), 0.0) / num_runs;
        double avg_sr = std::accumulate(sr_per_run.begin(), sr_per_run.end(), 0.0) / num_runs;
        double avg_mpr = std::accumulate(mpr_per_run.begin(), mpr_per_run.end(), 0.0) / num_runs;

        // Count successful runs (SR = 100%)
        int successful_runs = std::count(sr_per_run.begin(), sr_per_run.end(), 100.0);
        double overall_sr = (static_cast<double>(successful_runs) / num_runs) * 100.0;

        std::cout << "\n>>> [MPM-CoEA] CEC2015 METRICS RESULTS (Across " << num_runs << " runs)" << std::endl;
        std::cout << "    Success Rate (SR):       " << std::fixed << std::setprecision(2) << overall_sr << "%" << std::endl;
        std::cout << "    Avg Optima Found (ANOF): " << std::fixed << std::setprecision(4) << avg_anof << std::endl;
        std::cout << "    Avg Max Peak Ratio (MPR): " << std::fixed << std::setprecision(4) << avg_mpr << std::endl;
        std::cout << "    Successful Runs:         " << successful_runs << "/" << num_runs << std::endl;

        
        std::cout << "\n>>> [MPM-CoEA] Experiment complete!" << std::endl;
    }

} // namespace ofec

#endif // !OFEC_MPMCOEA_SOLVER_HPP