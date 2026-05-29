#ifndef OFEC_MPMCOEA_SOLVER_HPP
#define OFEC_MPMCOEA_SOLVER_HPP

#include "custom_method.hpp"
#include "../instance/problem/continuous/free_peaks/free_peaks.h"
#include "../instance/problem/continuous/free_peaks/subproblem/transform/transform_x/transform_x_base.h"
#include "../instance/problem/continuous/free_peaks/subproblem/transform/transform_y/transform_y_base.h"
#include "../instance/problem/continuous/free_peaks/factory.h"
#include "../instance/problem/continuous/free_peaks/subproblem/subproblem.h"
#include "../core/environment/environment.h"
#include <queue>
#include <numeric>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numbers>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <memory>
#include <random>
#include <filesystem>
#include <string>
#include <limits>
#include <cassert>
#include <sstream>
#include <functional>
#include <unordered_map>
#include <set>

namespace ofec {

    // NEAREST-BETTER CLUSTERING (NBC)
    std::vector<std::vector<int>> nearestBetterClustering(
        const std::vector<std::vector<double>>& population,
        const std::vector<double>& fitness,
        double phi = 2.0)
    {
        const int n = static_cast<int>(population.size());
        if (n == 0) return {};

        struct Edge { int from, to; double weight; };
        std::vector < Edge > edges;
        edges.reserve(static_cast<size_t>(n));

        for (int i = 0; i < n; i++) {
            double min_dist = std::numeric_limits<double>::max();
            int nb_idx = -1;
            for (int j = 0; j < n; j++) {
                if (fitness[j] > fitness[i] + 1e-12) {
                    double d = 0.0;
                    for (size_t k = 0; k < population[i].size(); k++) {
                        const double diff = population[i][k] - population[j][k];
                        d += diff * diff;
                    }
                    d = std::sqrt(d);
                    if (d < min_dist) { min_dist = d; nb_idx = j; }
                }
            }
            if (nb_idx != -1) edges.push_back({ i, nb_idx, min_dist });
        }

        if (edges.empty()) {
            std::vector<std::vector<int>> clusters;
            clusters.reserve(static_cast<size_t>(n));
            for (int i = 0; i < n; i++) clusters.push_back({ i });
            return clusters;
        }

        double sum_w = 0.0;
        for (const auto& e : edges) sum_w += e.weight;
        const double threshold = phi * (sum_w / static_cast<double>(edges.size()));

        std::vector<std::vector<int>> adj(static_cast<size_t>(n));
        for (const auto& e : edges)
            if (e.weight <= threshold) {
                adj[static_cast<size_t>(e.from)].push_back(e.to);
                adj[static_cast<size_t>(e.to)].push_back(e.from);
            }

        std::vector<std::vector<int>> clusters;
        std::vector<bool> visited(static_cast<size_t>(n), false);
        for (int i = 0; i < n; i++) {
            if (!visited[static_cast<size_t>(i)]) {
                std::vector<int> cluster;
                std::queue<int> q;
                q.push(i); visited[static_cast<size_t>(i)] = true;
                while (!q.empty()) {
                    int u = q.front(); q.pop();
                    cluster.push_back(u);
                    for (int v : adj[static_cast<size_t>(u)])
                        if (!visited[static_cast<size_t>(v)]) {
                            visited[static_cast<size_t>(v)] = true;
                            q.push(v);
                        }
                }
                clusters.push_back(std::move(cluster));
            }
        }
        return clusters;
    }


    // MPMMO_BENCHMARK
    class MPMMO_Benchmark {
    private:
        FreePeaks* m_problem;
        Environment* m_env;
        int m_num_dms;
        std::vector<std::vector<int>> m_party_indices;

    public:
        MPMMO_Benchmark(FreePeaks* prob, int num_dms, Environment* env)
            : m_problem(prob), m_env(env), m_num_dms(num_dms) {
            decomposeVariables();
        }

        void decomposeVariables() {
            const int total_vars = static_cast<int>(m_problem->numberVariables());
            m_party_indices.resize(static_cast<size_t>(m_num_dms));
            for (int v = 0; v < total_vars; v++)
                m_party_indices[static_cast<size_t>(v % m_num_dms)].push_back(v);

            std::cout << "[INFO] Variable decomposition (dim=" << total_vars
                << ", parties=" << m_num_dms << "):\n";
            for (int dm = 0; dm < m_num_dms; dm++) {
                std::cout << "  Party " << dm << " -> vars [";
                for (size_t i = 0; i < m_party_indices[dm].size(); i++) {
                    std::cout << m_party_indices[dm][i];
                    if (i + 1 < m_party_indices[dm].size()) std::cout << ", ";
                }
                std::cout << "]\n";
            }
        }

        double evaluateJoint(const std::vector<double>& full_sol, Environment* env) const {
            std::unique_ptr<SolutionBase> sol(m_problem->createSolution());
            auto& sv = dynamic_cast <Solution < VariableVector < Real > > &> (*sol).variable().vector();
            for (size_t i = 0; i < sv.size() && i < full_sol.size(); i++)
                sv[i] = static_cast<Real>(full_sol[i]);
            sol->evaluate(env, false);
            double result = static_cast<double>(sol->objective(0));
            return std::isfinite(result) ? result : -1e30;
        }

        double evaluateDM(int dm_id, const std::vector<double>& my_vars,
            const std::vector<std::vector<double>>& context, Environment* env) const {
            const int total_dim = static_cast<int>(m_problem->numberVariables());
            std::vector<double> full_sol(static_cast<size_t>(total_dim), 0.5);
            const auto& my_idx = m_party_indices[static_cast<size_t>(dm_id)];
            for (size_t i = 0; i < my_idx.size() && i < my_vars.size(); i++)
                full_sol[static_cast<size_t>(my_idx[i])] = my_vars[i];

            for (int other = 0; other < m_num_dms; other++) {
                if (other == dm_id) continue;
                const auto& oth_idx = m_party_indices[static_cast<size_t>(other)];
                const auto& ctx = context[static_cast<size_t>(other)];
                for (size_t i = 0; i < oth_idx.size() && i < ctx.size(); i++)
                    full_sol[static_cast<size_t>(oth_idx[i])] = ctx[i];
            }
            return evaluateJoint(full_sol, env);
        }

        int getNumDMs() const { return m_num_dms; }
        int getNumVariables() const { return static_cast<int>(m_problem->numberVariables()); }
        int getPartyDim(int p) const { return static_cast<int>(m_party_indices[static_cast<size_t>(p)].size()); }
        const std::vector<int>& getPartyIndices(int p) const { return m_party_indices[static_cast<size_t>(p)]; }
        FreePeaks* getProblem() const { return m_problem; }
        Environment* getEnv() const { return m_env; }
    };


    // MPM_COEA
    class MPM_CoEA {
    private:
        std::shared_ptr < MPMMO_Benchmark > m_bench;
        std::shared_ptr < Environment > m_env;

        // Algorithm parameters
        int m_pop_size;
        int m_med_pop_size;
        int m_K;

        static constexpr int MIN_POP_SIZE = 20;
        static constexpr int STAGNATION_LIMIT = 80;

        // ADAPTIVE DE PARAMETERS - vary by phase
        double m_de_F = 0.65;
        double m_de_CR = 0.85;

        // FITNESS SHARING PARAMETERS
        // Increased sharing radius for 4D space (was 0.15)
        static constexpr double SHARING_RADIUS = 0.35;
        static constexpr double SHARING_ALPHA = 1.0;

        // Competing populations
        std::vector<std::vector<std::vector<std::vector<double>>>> m_competing_pops;
        std::vector<std::vector<std::vector<double>>> m_competing_fitness;

        // Mediating population
        std::vector<std::vector<double>> m_mediating_pop;
        std::vector<double> m_mediating_fitness;

        // MULTI-CONTEXT ARCHIVE
        std::vector<std::vector<std::vector<double>>> m_context_archive;
        std::vector<std::vector<double>> m_context_fitness;

        // PER-SUB-POP CONTEXT ASSIGNMENT
        std::vector<std::vector<std::vector<int>>> m_subpop_context_assignment;

        // FOUND OPTIMA ARCHIVE
        struct FoundOptimum {
            std::vector<double> sol;
            double fitness;
            int basin_id = -1;
        };
        std::vector < FoundOptimum > m_found_optima;
        static constexpr double CLEARING_R = 0.20;

        // Added new Tabu peak-exclusion mechanism
        struct TabuEntry {
            std::vector<double> center;
            int expire_gen;
        };
        std::vector<TabuEntry> m_tabu_peaks;

        // Elite archive
        std::vector<std::vector<std::vector<double>>> m_elite_archive;
        std::vector<std::vector<double>> m_elite_fitness;

        // Evaluation caches
        mutable std::unordered_map<std::string, double> m_cache_joint;
        mutable std::unordered_map<std::string, double> m_cache_dm;

        static std::string vecKey(const std::vector<double>& v) {
            std::string s;
            s.reserve(v.size() * sizeof(double));
            for (double d : v) s.append(reinterpret_cast<const char*>(&d), sizeof(d));
            return s;
        }

        double evalJointCached(const std::vector<double>& sol) const {
            auto key = vecKey(sol);
            auto it = m_cache_joint.find(key);
            if (it != m_cache_joint.end()) return it->second;
            double f = m_bench->evaluateJoint(sol, m_env.get());
            m_cache_joint[key] = f;
            return f;
        }

        double evalDMCached(int dm, const std::vector<double>& my_vars,
            const std::vector<std::vector<double>>& ctx) const {
            std::string key;
            key.reserve(sizeof(int) + my_vars.size() * sizeof(double) + ctx.size() * ctx[0].size() * sizeof(double));
            key.append(reinterpret_cast<const char*>(&dm), sizeof(dm));
            key.append(vecKey(my_vars));
            for (const auto& c : ctx) key.append(vecKey(c));
            auto it = m_cache_dm.find(key);
            if (it != m_cache_dm.end()) return it->second;
            double f = m_bench->evaluateDM(dm, my_vars, ctx, m_env.get());
            m_cache_dm[key] = f;
            return f;
        }

        void clearCaches() const {
            if (m_cache_joint.size() > 200000) m_cache_joint.clear();
            if (m_cache_dm.size() > 200000) m_cache_dm.clear();
        }

        std::string m_output_directory;
        int m_generation = 0;
        int m_stagnation_counter = 0;
        double m_best_prev = -std::numeric_limits<double>::max();
        std::vector<std::pair<int, double>> m_convergence_log;

        // Phase tracking
        // Added new: Immediate cooperative phase (no neutral burn-in)
        bool m_neutral_phase = false;
        int m_neutral_gens = 0;

        // Restart tracking for adaptive strategy
        int m_restart_count = 0;
        std::vector<double> m_restart_best_history;

    public:
        MPM_CoEA(std::shared_ptr < MPMMO_Benchmark > bench,
            std::shared_ptr < Environment > env,
            int pop_size = 40,
            int med_pop = 400,
            int K = 10,
            const std::string& output_dir = "visualization/solutions")
            : m_bench(bench), m_env(env),
            m_pop_size(pop_size), m_med_pop_size(med_pop), m_K(K),
            m_output_directory(output_dir) {
            std::filesystem::create_directories(output_dir);
        }

        // Added new Tabu helpers
        bool isTabu(const std::vector<double>& sol) const {
            for (const auto& t : m_tabu_peaks) {
                if (m_generation > t.expire_gen) continue;
                double d2 = 0.0;
                for (size_t d = 0; d < sol.size(); ++d) {
                    double diff = sol[d] - t.center[d];
                    d2 += diff * diff;
                }
                if (std::sqrt(d2) < 0.30) return true;
            }
            return false;
        }

        void addTabuPeak(const std::vector<double>& sol) {
            m_tabu_peaks.push_back({ sol, m_generation + 80 });
        }

        // Added new: Adaptive clearing radius from peak statistics
        double getAdaptiveClearingRadius() const {
            if (m_found_optima.size() < 2) return 0.30;
            double min_dist = 1e30;
            for (size_t i = 0; i < m_found_optima.size(); ++i) {
                for (size_t j = i + 1; j < m_found_optima.size(); ++j) {
                    double d2 = 0.0;
                    for (size_t d = 0; d < m_found_optima[i].sol.size(); ++d) {
                        double diff = m_found_optima[i].sol[d] - m_found_optima[j].sol[d];
                        d2 += diff * diff;
                    }
                    min_dist = std::min(min_dist, std::sqrt(d2));
                }
            }
            double r = min_dist * 0.45;
            if (r < 0.20) r = 0.20;
            if (r > 0.50) r = 0.50;
            return r;
        }

        // INITIALIZATION
        void initialize() {
            m_cache_joint.clear();
            m_cache_dm.clear();
            m_tabu_peaks.clear();  // Added new: clear tabu at start
            const int total_dim = m_bench->getNumVariables();
            const int num_dms = m_bench->getNumDMs();

            // Stratified pool with Latin Hypercube Sampling for better coverage
            const int POOL_SIZE = 16000;

            std::mt19937 rng(42u);
            std::uniform_real_distribution<double> unif(0.0, 1.0);

            // Latin Hypercube Sampling for initial pool
            std::vector<std::vector<double>> pool(
                static_cast<size_t>(POOL_SIZE),
                std::vector<double>(static_cast<size_t>(total_dim)));

            for (int d = 0; d < total_dim; ++d) {
                std::vector<double> strata(POOL_SIZE);
                for (int i = 0; i < POOL_SIZE; ++i) {
                    strata[i] = (i + unif(rng)) / POOL_SIZE;
                }
                std::shuffle(strata.begin(), strata.end(), rng);
                for (int i = 0; i < POOL_SIZE; ++i) {
                    pool[i][d] = strata[i];
                }
            }

            // Joint fitness evaluation
            std::vector<double> joint_fitness(static_cast<size_t>(POOL_SIZE));
            for (int i = 0; i < POOL_SIZE; i++)
                joint_fitness[static_cast<size_t>(i)] = evalJointCached(pool[static_cast<size_t>(i)]);

            {
                double mx = *std::max_element(joint_fitness.begin(), joint_fitness.end());
                double mn = *std::min_element(joint_fitness.begin(), joint_fitness.end());
                std::cout << "[INFO] Pool joint fitness range: ["
                    << std::fixed << std::setprecision(3)
                    << mn << ", " << mx << "]\n";
            }

            // Neutral context
            std::vector<std::vector<double>> neutral_ctx(static_cast<size_t>(num_dms));
            for (int dm = 0; dm < num_dms; dm++)
                neutral_ctx[static_cast<size_t>(dm)].assign(
                    static_cast<size_t>(m_bench->getPartyDim(dm)), 0.5);

            // Per-party fitness under neutral context
            std::vector<std::vector<double>> per_fit(
                static_cast<size_t>(num_dms),
                std::vector<double>(static_cast<size_t>(POOL_SIZE)));

            for (int dm = 0; dm < num_dms; dm++) {
                const int dm_dim = m_bench->getPartyDim(dm);
                const auto& idx = m_bench->getPartyIndices(dm);
                for (int i = 0; i < POOL_SIZE; i++) {
                    std::vector<double> my_vars(static_cast<size_t>(dm_dim));
                    for (int d = 0; d < dm_dim; d++)
                        my_vars[static_cast<size_t>(d)] = pool[static_cast<size_t>(i)][static_cast<size_t>(idx[d])];
                    per_fit[static_cast<size_t>(dm)][static_cast<size_t>(i)] =
                        evalDMCached(dm, my_vars, neutral_ctx);
                }
            }

            // NBC clustering per party
            m_competing_pops.resize(static_cast<size_t>(num_dms));
            m_competing_fitness.resize(static_cast<size_t>(num_dms));
            m_elite_archive.resize(static_cast<size_t>(num_dms));
            m_elite_fitness.resize(static_cast<size_t>(num_dms));

            for (int dm = 0; dm < num_dms; dm++) {
                const int dm_dim = m_bench->getPartyDim(dm);
                const auto& idx = m_bench->getPartyIndices(dm);

                std::vector<std::vector<double>> dm_pool(
                    static_cast<size_t>(POOL_SIZE),
                    std::vector<double>(static_cast<size_t>(dm_dim)));
                for (int i = 0; i < POOL_SIZE; i++)
                    for (int d = 0; d < dm_dim; d++)
                        dm_pool[static_cast<size_t>(i)][static_cast<size_t>(d)] =
                        pool[static_cast<size_t>(i)][static_cast<size_t>(idx[d])];

                auto clusters = nearestBetterClustering(dm_pool, per_fit[static_cast<size_t>(dm)], 2.0);

                std::sort(clusters.begin(), clusters.end(),
                    [&](const std::vector<int>& a, const std::vector<int>& b) {
                        double ba = -1e30, bb = -1e30;
                        for (int ii : a) ba = std::max(ba, per_fit[dm][static_cast<size_t>(ii)]);
                        for (int ii : b) bb = std::max(bb, per_fit[dm][static_cast<size_t>(ii)]);
                        return ba > bb;
                    });

                if (static_cast<int>(clusters.size()) > m_K)
                    clusters.resize(static_cast<size_t>(m_K));

                m_competing_pops[static_cast<size_t>(dm)].clear();
                m_competing_fitness[static_cast<size_t>(dm)].clear();

                std::uniform_real_distribution<double> noise(-0.02, 0.02);

                for (auto& cluster : clusters) {
                    std::sort(cluster.begin(), cluster.end(),
                        [&](int a, int b) {
                            return per_fit[dm][static_cast<size_t>(a)] > per_fit[dm][static_cast<size_t>(b)];
                        });

                    int take = std::min(static_cast<int>(cluster.size()), m_pop_size);

                    std::vector<std::vector<double>> sub_pop(
                        static_cast<size_t>(take),
                        std::vector<double>(static_cast<size_t>(dm_dim)));
                    std::vector<double> sub_fit(static_cast<size_t>(take));

                    for (int k = 0; k < take; k++) {
                        const int pi = cluster[static_cast<size_t>(k)];
                        for (int d = 0; d < dm_dim; d++)
                            sub_pop[static_cast<size_t>(k)][static_cast<size_t>(d)] =
                            dm_pool[static_cast<size_t>(pi)][static_cast<size_t>(d)];
                        sub_fit[static_cast<size_t>(k)] = per_fit[dm][static_cast<size_t>(pi)];
                    }

                    while (static_cast<int>(sub_pop.size()) < MIN_POP_SIZE) {
                        std::vector<double> perturbed = sub_pop[0];
                        for (auto& v : perturbed) {
                            v += noise(rng);
                            v = std::max(0.0, std::min(1.0, v));
                        }
                        sub_pop.push_back(perturbed);
                        sub_fit.push_back(
                            m_bench->evaluateDM(dm, perturbed, neutral_ctx, m_env.get()));
                    }

                    m_competing_pops[static_cast<size_t>(dm)].push_back(std::move(sub_pop));
                    m_competing_fitness[static_cast<size_t>(dm)].push_back(std::move(sub_fit));
                }

                buildEliteArchive(dm, dm_pool, per_fit[dm]);

                std::cout << "[INFO] Party " << dm << ": "
                    << m_competing_pops[dm].size() << " sub-pop(s), sizes:";
                for (const auto& sp : m_competing_pops[dm])
                    std::cout << " " << sp.size();
                std::cout << "\n";
            }

            // Mediating population - diverse selection
            std::vector<int> joint_order(static_cast<size_t>(POOL_SIZE));
            std::iota(joint_order.begin(), joint_order.end(), 0);
            std::sort(joint_order.begin(), joint_order.end(),
                [&](int a, int b) {
                    return joint_fitness[static_cast<size_t>(a)] > joint_fitness[static_cast<size_t>(b)];
                });

            m_mediating_pop.clear();
            m_mediating_fitness.clear();
            m_mediating_pop.reserve(static_cast<size_t>(m_med_pop_size));
            m_mediating_fitness.reserve(static_cast<size_t>(m_med_pop_size));

            const double diversity_r = 0.08;
            for (int rank_i : joint_order) {
                if (static_cast<int>(m_mediating_pop.size()) >= m_med_pop_size / 2) break;
                const auto& cand = pool[static_cast<size_t>(rank_i)];
                bool too_close = false;
                for (const auto& existing : m_mediating_pop) {
                    double d2 = 0.0;
                    for (size_t dd = 0; dd < cand.size(); dd++) {
                        double diff = cand[dd] - existing[dd];
                        d2 += diff * diff;
                    }
                    if (std::sqrt(d2) < diversity_r) { too_close = true; break; }
                }
                if (!too_close) {
                    m_mediating_pop.push_back(cand);
                    m_mediating_fitness.push_back(joint_fitness[static_cast<size_t>(rank_i)]);
                }
            }

            while (static_cast<int>(m_mediating_pop.size()) < m_med_pop_size) {
                std::vector<double> rand_sol(static_cast<size_t>(total_dim));
                for (int d = 0; d < total_dim; ++d) rand_sol[d] = unif(rng);
                m_mediating_pop.push_back(rand_sol);
                m_mediating_fitness.push_back(evalJointCached(rand_sol));
            }

            // Initialize context archives with diversity enforcement
            const int num_dms_init = m_bench->getNumDMs();
            m_context_archive.resize(num_dms_init);
            m_context_fitness.resize(num_dms_init);

            // Initialize per-sub-pop context assignment matrix
            m_subpop_context_assignment.resize(num_dms_init);

            for (int dm = 0; dm < num_dms_init; ++dm) {
                m_context_archive[dm].clear();
                m_context_fitness[dm].clear();
                for (size_t sp = 0; sp < m_competing_pops[dm].size(); ++sp) {
                    if (m_competing_pops[dm][sp].empty()) continue;
                    m_context_archive[dm].push_back(m_competing_pops[dm][sp][0]);
                    m_context_fitness[dm].push_back(m_competing_fitness[dm][sp][0]);
                }
                pruneContextArchiveDiverse(dm, 0.20);

                // Initialize context assignment for this party's sub-pops
                const int n_subpops = static_cast<int>(m_competing_pops[dm].size());
                m_subpop_context_assignment[dm].resize(n_subpops);
                for (int sp = 0; sp < n_subpops; ++sp) {
                    m_subpop_context_assignment[dm][sp].resize(num_dms_init, -1);
                }
            }

            // Initialize found optima
            updateFoundOptimaArchive();
            pruneFoundOptimaArchive();

            m_best_prev = getBestConsensusFitness();
            m_stagnation_counter = 0;
            m_convergence_log.clear();
            m_convergence_log.push_back({ 0, m_best_prev });
            m_restart_count = 0;

            std::cout << "[INFO] Init complete | med_pop=" << m_mediating_pop.size()
                << " | best_joint=" << std::fixed << std::setprecision(4)
                << m_best_prev << "\n";
        }


        // RUN ONE GENERATION
        void run_one_generation() {
            ++m_generation;
            clearCaches();
            const int num_dms = m_bench->getNumDMs();
            const int total_dim = m_bench->getNumVariables();

            // Added new: No neutral phase; always use cooperative parameters
            m_de_F = 0.55;
            m_de_CR = 0.95;

            // DE on competing sub-populations with DISTINCT contexts
            for (int dm = 0; dm < num_dms; dm++) {
                const int dm_dim = m_bench->getPartyDim(dm);

                for (size_t sp = 0; sp < m_competing_pops[dm].size(); sp++) {
                    auto& pop = m_competing_pops[dm][sp];
                    auto& fit = m_competing_fitness[dm][sp];
                    const int pop_n = static_cast<int>(pop.size());
                    if (pop_n < 4) continue;

                    // Build per-sub-pop context with DIVERSITY CONSTRAINT
                    std::vector<std::vector<double>> ctx(num_dms);
                    for (int other = 0; other < num_dms; ++other) {
                        if (other == dm) {
                            ctx[other].assign(dm_dim, 0.5);
                        }
                        else {
                            const auto& arch = m_context_archive[other];
                            if (!arch.empty()) {
                                int ctx_idx = getDiverseContextIndex(dm, static_cast<int>(sp), other);
                                ctx[other] = arch[ctx_idx % arch.size()];
                                m_subpop_context_assignment[dm][sp][other] = ctx_idx;
                            }
                            else {
                                ctx[other].assign(m_bench->getPartyDim(other), 0.5);
                            }
                        }
                    }

                    std::mt19937 rng_de(static_cast<unsigned>(
                        m_generation * 100003u + dm * 1009u + static_cast<unsigned>(sp)));
                    std::uniform_real_distribution<double> unif(0.0, 1.0);
                    std::uniform_int_distribution<int> pick(0, pop_n - 1);

                    // Added new: Generate all trials first, then tournament-select survivors
                    std::vector<std::vector<double>> trials(pop_n);
                    std::vector<double> trial_fits(pop_n, -1e30);
                    std::vector<bool> trial_valid(pop_n, false);

                    for (int i = 0; i < pop_n; i++) {
                        int r1 = i, r2 = i, r3 = i, att = 0;
                        while (r1 == i && ++att < 2000) r1 = pick(rng_de);
                        att = 0; while ((r2 == i || r2 == r1) && ++att < 2000) r2 = pick(rng_de);
                        att = 0; while ((r3 == i || r3 == r1 || r3 == r2) && ++att < 2000) r3 = pick(rng_de);
                        if (r1 == i || r2 == i || r3 == i) continue;

                        std::vector<double> trial(pop[i]);
                        const int jrand = static_cast<int>(unif(rng_de) * dm_dim);
                        for (int d = 0; d < dm_dim; d++) {
                            if (unif(rng_de) < m_de_CR || d == jrand) {
                                double v = pop[r1][d] + m_de_F * (pop[r2][d] - pop[r3][d]);
                                v = std::max(0.0, std::min(1.0, v));
                                trial[d] = v;
                            }
                        }

                        const double tf = evalDMCached(dm, trial, ctx);
                        if (std::isfinite(tf)) {
                            trials[i] = std::move(trial);
                            trial_fits[i] = tf;
                            trial_valid[i] = true;
                        }
                    }

                    // Added new: Tournament-size-2 survival from merged parent+offspring pool
                    std::vector<std::vector<double>> new_pop;
                    std::vector<double> new_fit;
                    new_pop.reserve(pop_n);
                    new_fit.reserve(pop_n);

                    for (int i = 0; i < pop_n; ++i) {
                        // Parent candidate
                        std::vector<double> cand_p = pop[i];
                        double fit_p = fit[i];
                        // Offspring candidate
                        std::vector<double> cand_o = trial_valid[i] ? trials[i] : pop[i];
                        double fit_o = trial_valid[i] ? trial_fits[i] : fit[i];

                        // Added new: Tabu check - if offspring is tabu, force parent
                        if (trial_valid[i] && isTabu(cand_o)) {
                            cand_o = pop[i];
                            fit_o = fit[i];
                        }

                        // Tournament size 2
                        std::uniform_int_distribution<int> coin(0, 1);
                        bool use_offspring = (fit_o > fit_p);
                        if (std::abs(fit_o - fit_p) < 1e-12) {
                            use_offspring = (coin(rng_de) == 1);
                        }

                        if (use_offspring) {
                            new_pop.push_back(std::move(cand_o));
                            new_fit.push_back(fit_o);
                        }
                        else {
                            new_pop.push_back(std::move(cand_p));
                            new_fit.push_back(fit_p);
                        }
                    }

                    pop = std::move(new_pop);
                    fit = std::move(new_fit);
                }
            }

            // UPDATE CONTEXT ARCHIVES with diversity enforcement
            for (int dm = 0; dm < num_dms; dm++) {
                std::vector<std::vector<double>> candidates;
                std::vector<double> cand_fit;

                // Sub-population bests
                for (size_t sp = 0; sp < m_competing_pops[dm].size(); sp++) {
                    if (m_competing_pops[dm][sp].empty()) continue;
                    candidates.push_back(m_competing_pops[dm][sp][0]);
                    cand_fit.push_back(m_competing_fitness[dm][sp][0]);
                }

                // Existing archive
                for (size_t i = 0; i < m_context_archive[dm].size(); ++i) {
                    candidates.push_back(m_context_archive[dm][i]);
                    cand_fit.push_back(m_context_fitness[dm][i]);
                }

                // Top mediating projections
                std::vector<size_t> order(m_mediating_pop.size());
                std::iota(order.begin(), order.end(), 0);
                std::sort(order.begin(), order.end(),
                    [&](size_t a, size_t b) { return m_mediating_fitness[a] > m_mediating_fitness[b]; });
                const auto& idx = m_bench->getPartyIndices(dm);
                const int dm_dim = m_bench->getPartyDim(dm);
                for (size_t r = 0; r < order.size() && r < static_cast<size_t>(m_K); ++r) {
                    std::vector<double> ctx(dm_dim);
                    for (int d = 0; d < dm_dim; ++d) ctx[d] = m_mediating_pop[order[r]][idx[d]];
                    candidates.push_back(std::move(ctx));
                    cand_fit.push_back(m_mediating_fitness[order[r]]);
                }

                // Random mediating projections
                std::mt19937 rng_div(static_cast<unsigned>(m_generation * 44444u + dm));
                std::uniform_int_distribution<int> pick_rand(0, static_cast<int>(m_mediating_pop.size()) - 1);
                for (int rand_i = 0; rand_i < 20; ++rand_i) {
                    int ri = pick_rand(rng_div);
                    std::vector<double> ctx(dm_dim);
                    for (int d = 0; d < dm_dim; ++d) ctx[d] = m_mediating_pop[ri][idx[d]];
                    candidates.push_back(std::move(ctx));
                    cand_fit.push_back(m_mediating_fitness[ri]);
                }

                m_context_archive[dm] = std::move(candidates);
                m_context_fitness[dm] = std::move(cand_fit);

                // STRONGER diversity pruning
                pruneContextArchiveDiverse(dm, 0.15);

                updateEliteArchive(dm);
            }

            // SPECIES-BASED MEDIATING EVOLUTION with FITNESS SHARING
            evolveMediatingSpeciesFitnessSharing();

            // Update explicit optima archive
            updateFoundOptimaArchive();
            pruneFoundOptimaArchive();

            if (m_generation % 20 == 0) {
                std::cout << "  [archive fitnesses:";
                for (const auto& opt : m_found_optima) {
                    std::cout << " " << std::fixed << std::setprecision(1) << opt.fitness;
                }
                std::cout << "]\n";
            }

            // Inject consensus candidates
            injectConsensusCandidates();

            // Enforce per-basin capacity
            clearMediatingPopulation();

            // Feed diverse mediating projections back
            mergeContextsFromMediating();

            // Periodic cooperative exploration
            if (m_generation % 5 == 0) {
                injectCartesianFromSubPops();
            }

            // STAGNATION DETECTION with ADAPTIVE RESTART
            const double cur_best = getBestConsensusFitness();
            m_convergence_log.push_back({ m_generation, cur_best });

            if (cur_best > m_best_prev + 1e-6) {
                m_best_prev = cur_best;
                m_stagnation_counter = 0;
            }
            else {
                ++m_stagnation_counter;
            }

            if (m_stagnation_counter >= STAGNATION_LIMIT) {
                std::cout << "[RESTART] Stagnation at gen " << m_generation
                    << " (best=" << std::fixed << std::setprecision(4)
                    << cur_best << ") — FULL RESET";

                ++m_restart_count;
                m_restart_best_history.push_back(cur_best);

                // ADAPTIVE RESTART: Full context archive reset
                adaptiveRestart();
                m_stagnation_counter = 0;
            }

            if (m_generation % 20 == 0) {
                std::cout << "[gen " << std::setw(4) << m_generation
                    << "] best_joint=" << std::fixed << std::setprecision(4)
                    << cur_best << " | stagnation=" << m_stagnation_counter
                    << " | found=" << m_found_optima.size()
                    << " | arch:";
                for (int dm = 0; dm < num_dms; dm++) {
                    std::cout << " DM" << dm << "=" << m_context_archive[dm].size();
                }
                std::cout << "\n";
            }
        }


        // DISTINCT CONTEXT INDEX - Ensures diversity
        int getDiverseContextIndex(int dm, int sp, int other_dm) {
            const auto& arch = m_context_archive[other_dm];
            if (arch.empty()) return 0;

            std::set<int> used_indices;
            for (size_t other_sp = 0; other_sp < m_subpop_context_assignment[dm].size(); ++other_sp) {
                if (static_cast<int>(other_sp) != sp &&
                    other_dm < static_cast<int>(m_subpop_context_assignment[dm][other_sp].size()) &&
                    m_subpop_context_assignment[dm][other_sp][other_dm] >= 0) {
                    used_indices.insert(m_subpop_context_assignment[dm][other_sp][other_dm] % static_cast<int>(arch.size()));
                }
            }

            for (int idx = 0; idx < static_cast<int>(arch.size()) * 2; ++idx) {
                int actual_idx = idx % static_cast<int>(arch.size());
                if (used_indices.find(actual_idx) == used_indices.end()) {
                    return actual_idx;
                }
            }

            return sp % static_cast<int>(arch.size());
        }

        void rebuildContextAssignments() {
            const int num_dms = m_bench->getNumDMs();
            for (int dm = 0; dm < num_dms; ++dm) {
                const int n_subpops = static_cast<int>(m_competing_pops[dm].size());
                m_subpop_context_assignment[dm].resize(n_subpops);
                for (int sp = 0; sp < n_subpops; ++sp) {
                    m_subpop_context_assignment[dm][sp].assign(num_dms, -1);
                }
            }
        }


        // FITNESS SHARING in Mediating Population
        double fitnessSharingPenalty(const std::vector<double>& sol, double raw_fitness) const {
            if (m_found_optima.empty()) return raw_fitness;

            double niche_count = 1.0;
            for (const auto& opt : m_found_optima) {
                if (opt.fitness < 30.0) continue;

                double d2 = 0.0;
                for (size_t d = 0; d < sol.size(); ++d) {
                    double diff = sol[d] - opt.sol[d];
                    d2 += diff * diff;
                }
                double dist = std::sqrt(d2);

                if (dist < SHARING_RADIUS) {
                    double sh = 1.0 - std::pow(dist / SHARING_RADIUS, SHARING_ALPHA);
                    niche_count += std::max(0.0, sh);
                }
            }

            return raw_fitness / niche_count;
        }

        // EVOLVE MEDIATING SPECIES with FITNESS SHARING
        void evolveMediatingSpeciesFitnessSharing() {
            const int total_dim = m_bench->getNumVariables();
            const int med_n = static_cast<int>(m_mediating_pop.size());
            if (med_n < 8) return;

            // Sample for NBC
            std::vector<std::vector<double>> sample_pop;
            std::vector<double> sample_fit;
            if (med_n > 300) {
                std::vector<size_t> order(med_n);
                std::iota(order.begin(), order.end(), 0);
                std::sort(order.begin(), order.end(),
                    [&](size_t a, size_t b) { return m_mediating_fitness[a] > m_mediating_fitness[b]; });
                for (int i = 0; i < 150; ++i) {
                    sample_pop.push_back(m_mediating_pop[order[i]]);
                    sample_fit.push_back(m_mediating_fitness[order[i]]);
                }
                std::mt19937 rng_samp(static_cast<unsigned>(m_generation * 33333u));
                std::uniform_int_distribution<int> pick(0, med_n - 1);
                for (int i = 0; i < 150; ++i) {
                    int idx = pick(rng_samp);
                    sample_pop.push_back(m_mediating_pop[idx]);
                    sample_fit.push_back(m_mediating_fitness[idx]);
                }
            }
            else {
                sample_pop = m_mediating_pop;
                sample_fit = m_mediating_fitness;
            }

            auto clusters = nearestBetterClustering(sample_pop, sample_fit, 1.5);

            std::vector<std::vector<int>> species;
            for (auto& c : clusters) {
                if (static_cast<int>(c.size()) >= 4) species.push_back(c);
            }
            while (species.size() < 2 && !species.empty()) {
                auto& largest = species[0];
                size_t half = largest.size() / 2;
                std::vector<int> second_half(largest.begin() + static_cast<int>(half), largest.end());
                largest.erase(largest.begin() + static_cast<int>(half), largest.end());
                species.push_back(second_half);
            }
            if (species.empty()) return;

            std::vector<size_t> sample_to_orig;
            if (med_n > 300) {
                std::vector<size_t> order(med_n);
                std::iota(order.begin(), order.end(), 0);
                std::sort(order.begin(), order.end(),
                    [&](size_t a, size_t b) { return m_mediating_fitness[a] > m_mediating_fitness[b]; });
                for (int i = 0; i < 150; ++i) sample_to_orig.push_back(order[i]);
                std::mt19937 rng_samp(static_cast<unsigned>(m_generation * 33333u));
                std::uniform_int_distribution<int> pick(0, med_n - 1);
                for (int i = 0; i < 150; ++i) sample_to_orig.push_back(pick(rng_samp));
            }
            else {
                for (size_t i = 0; i < static_cast<size_t>(med_n); ++i) sample_to_orig.push_back(i);
            }

            std::mt19937 rng_spec(static_cast<unsigned>(m_generation * 777777u));
            std::uniform_real_distribution<double> unif(0.0, 1.0);

            for (const auto& sp : species) {
                const int sp_size = static_cast<int>(sp.size());
                if (sp_size < 4) continue;

                std::uniform_int_distribution<int> pick_sp(0, sp_size - 1);

                for (int idx = 0; idx < sp_size; ++idx) {
                    int sample_i = sp[idx];
                    if (sample_i < 0 || sample_i >= static_cast<int>(sample_to_orig.size())) continue;
                    size_t global_i = sample_to_orig[static_cast<size_t>(sample_i)];
                    int r1 = idx, r2 = idx, r3 = idx, att = 0;
                    while (r1 == idx && ++att < 2000) r1 = pick_sp(rng_spec);
                    att = 0; while ((r2 == idx || r2 == r1) && ++att < 2000) r2 = pick_sp(rng_spec);
                    att = 0; while ((r3 == idx || r3 == r1 || r3 == r2) && ++att < 2000) r3 = pick_sp(rng_spec);
                    if (r1 == idx || r2 == idx || r3 == idx) continue;

                    size_t g1 = sample_to_orig[static_cast<size_t>(sp[r1])];
                    size_t g2 = sample_to_orig[static_cast<size_t>(sp[r2])];
                    size_t g3 = sample_to_orig[static_cast<size_t>(sp[r3])];

                    std::vector<double> trial(m_mediating_pop[global_i]);
                    const int jrand = static_cast<int>(unif(rng_spec) * total_dim);
                    for (int d = 0; d < total_dim; d++) {
                        if (unif(rng_spec) < m_de_CR || d == jrand) {
                            double v = m_mediating_pop[g1][d] + m_de_F * (m_mediating_pop[g2][d] - m_mediating_pop[g3][d]);
                            v = std::max(0.0, std::min(1.0, v));
                            trial[d] = v;
                        }
                    }

                    // Large-jump mutation
                    if (unif(rng_spec) < 0.20) {
                        std::normal_distribution<double> gauss(0.0, 0.08);
                        for (int d = 0; d < total_dim; ++d) {
                            trial[d] += gauss(rng_spec);
                            trial[d] = std::max(0.0, std::min(1.0, trial[d]));
                        }
                    }

                    const double raw_tf = evalJointCached(trial);
                    if (!std::isfinite(raw_tf)) continue;

                    // Added new: Reject tabu trials before sharing comparison
                    if (isTabu(trial)) continue;

                    const double shared_trial = fitnessSharingPenalty(trial, raw_tf);
                    const double shared_curr = fitnessSharingPenalty(m_mediating_pop[global_i], m_mediating_fitness[global_i]);

                    if (shared_trial >= shared_curr) {
                        m_mediating_pop[global_i] = std::move(trial);
                        m_mediating_fitness[global_i] = raw_tf;
                    }
                }
            }
        }


        // ADAPTIVE RESTART with full context reset
        void adaptiveRestart() {
            clearCaches();

            // Clear low-quality archive entries
            std::vector < FoundOptimum > kept;
            for (const auto& opt : m_found_optima) {
                if (opt.fitness >= 60.0) kept.push_back(opt);
            }
            m_found_optima = std::move(kept);

            // Added new: Rebuild tabu list from all retained optima so restarts avoid known basins
            m_tabu_peaks.clear();
            for (const auto& opt : m_found_optima) {
                addTabuPeak(opt.sol);
            }

            const int num_dms = m_bench->getNumDMs();
            std::mt19937 rng_rs(static_cast<unsigned>(m_generation * 55555u + m_restart_count * 9999u));
            std::uniform_real_distribution<double> unif(0.0, 1.0);

            // FULL CONTEXT ARCHIVE RESET
            for (int dm = 0; dm < num_dms; dm++) {
                m_context_archive[dm].clear();
                m_context_fitness[dm].clear();

                for (size_t sp = 0; sp < m_competing_pops[dm].size(); ++sp) {
                    if (m_competing_pops[dm][sp].empty()) continue;
                    m_context_archive[dm].push_back(m_competing_pops[dm][sp][0]);
                    m_context_fitness[dm].push_back(m_competing_fitness[dm][sp][0]);
                }

                pruneContextArchiveDiverse(dm, 0.20);
            }

            rebuildContextAssignments();

            for (int dm = 0; dm < num_dms; dm++) {
                const int dm_dim = m_bench->getPartyDim(dm);
                const auto& idx = m_bench->getPartyIndices(dm);

                for (size_t sp = 0; sp < m_competing_pops[dm].size(); ++sp) {
                    auto& pop = m_competing_pops[dm][sp];
                    auto& fit = m_competing_fitness[dm][sp];
                    if (pop.empty()) continue;

                    std::vector<double> best_saved = pop[0];
                    double best_fit_saved = fit[0];

                    for (size_t k = 1; k < pop.size(); ++k) {
                        bool valid = false;
                        int attempts = 0;
                        while (!valid && attempts < 200) {
                            for (auto& v : pop[k]) v = unif(rng_rs);
                            valid = true;
                            // Added new: Also avoid tabu basins during restart reinit
                            if (isTabu(pop[k])) valid = false;
                            for (const auto& opt : m_found_optima) {
                                double d2 = 0.0;
                                for (int d = 0; d < dm_dim; ++d) {
                                    double diff = pop[k][d] - opt.sol[idx[d]];
                                    d2 += diff * diff;
                                }
                                if (std::sqrt(d2) < CLEARING_R) { valid = false; break; }
                            }
                            ++attempts;
                        }
                        std::vector<std::vector<double>> neutral_ctx(num_dms);
                        for (int d = 0; d < num_dms; ++d)
                            neutral_ctx[d].assign(m_bench->getPartyDim(d), 0.5);
                        fit[k] = evalDMCached(dm, pop[k], neutral_ctx);
                    }
                    pop[0] = best_saved;
                    fit[0] = best_fit_saved;
                }

                while (static_cast<int>(m_competing_pops[dm].size()) < m_K) {
                    std::vector<std::vector<double>> sub_pop(MIN_POP_SIZE,
                        std::vector<double>(dm_dim));
                    std::vector<double> sub_fit(MIN_POP_SIZE);

                    for (int k = 0; k < MIN_POP_SIZE; ++k) {
                        bool valid = false;
                        int attempts = 0;
                        while (!valid && attempts < 200) {
                            for (auto& v : sub_pop[k]) v = unif(rng_rs);
                            valid = true;
                            // Added new: Tabu check for new sub-populations
                            if (isTabu(sub_pop[k])) valid = false;
                            for (const auto& opt : m_found_optima) {
                                if (opt.fitness < 60.0) continue;
                                double d2 = 0.0;
                                for (int d = 0; d < dm_dim; ++d) {
                                    double diff = sub_pop[k][d] - opt.sol[idx[d]];
                                    d2 += diff * diff;
                                }
                                if (std::sqrt(d2) < CLEARING_R) {
                                    valid = false; break;
                                }
                            }
                            ++attempts;
                        }

                        std::vector<std::vector<double>> diverse_ctx(num_dms);
                        std::vector<std::vector<double>> neutral_ctx(num_dms);
                        for (int d = 0; d < num_dms; ++d) {
                            neutral_ctx[d].assign(m_bench->getPartyDim(d), 0.5);
                            if (d == dm) {
                                diverse_ctx[d].assign(m_bench->getPartyDim(d), 0.5);
                            }
                            else {
                                if (!m_context_archive[d].empty()) {
                                    std::uniform_int_distribution<int> pick_ctx(
                                        0, static_cast<int>(m_context_archive[d].size()) - 1);
                                    diverse_ctx[d] = m_context_archive[d][pick_ctx(rng_rs)];
                                }
                                else {
                                    diverse_ctx[d].assign(m_bench->getPartyDim(d), 0.5);
                                }
                            }
                        }
                        double f_div = evalDMCached(dm, sub_pop[k], diverse_ctx);
                        double f_neu = evalDMCached(dm, sub_pop[k], neutral_ctx);
                        sub_fit[k] = std::max(f_div, f_neu);
                    }
                    m_competing_pops[dm].push_back(std::move(sub_pop));
                    m_competing_fitness[dm].push_back(std::move(sub_fit));
                }
            }

            // Inject random explorers into mediating population
            {
                const int total_dim = m_bench->getNumVariables();
                std::mt19937 rng_div(static_cast<unsigned>(m_generation * 99999u));
                std::uniform_real_distribution<double> unif_explore(0.0, 1.0);
                const int n_replace = static_cast<int>(m_mediating_pop.size() / 4);
                for (int r = 0; r < n_replace; ++r) {
                    std::vector<double> rsol(static_cast<size_t>(total_dim));
                    bool valid = false;
                    int attempts = 0;
                    while (!valid && attempts < 200) {
                        for (auto& v : rsol) v = unif_explore(rng_div);
                        valid = true;
                        // Added new: Tabu check for mediating explorers
                        if (isTabu(rsol)) valid = false;
                        for (const auto& opt : m_found_optima) {
                            double d2 = 0.0;
                            for (int d = 0; d < total_dim; ++d) {
                                double diff = rsol[d] - opt.sol[d];
                                d2 += diff * diff;
                            }
                            if (std::sqrt(d2) < 0.15) { valid = false; break; }
                        }
                        ++attempts;
                    }
                    double rf = m_bench->evaluateJoint(rsol, m_env.get());
                    forceInsertMediating(rsol, rf);
                }
            }
        }


        // STRICTER LOCAL OPTIMALITY CHECK
        bool strictLocalOptimalityCheck(const std::vector<double>& sol, double fitness) const {
            std::mt19937 rng_check(static_cast<unsigned>(m_generation * 77777u));
            std::normal_distribution<double> noise(0.0, 0.005);

            for (int check = 0; check < 50; ++check) {
                auto perturbed = sol;
                for (auto& v : perturbed) {
                    v += noise(rng_check);
                    v = std::max(0.0, std::min(1.0, v));
                }
                double pf = evalJointCached(perturbed);
                if (pf > fitness + 0.5) {
                    return false;
                }
            }
            return true;
        }


        // UPDATE FOUND OPTIMA ARCHIVE
        void updateFoundOptimaArchive() {
            for (size_t i = 0; i < m_mediating_pop.size(); ++i) {
                if (!std::isfinite(m_mediating_fitness[i])) continue;
                if (m_mediating_fitness[i] < 55.0) continue;

                // Added new: Use adaptive clearing radius
                double clear_r = getAdaptiveClearingRadius();

                double min_dist = 1e30;
                for (const auto& opt : m_found_optima) {
                    double d2 = 0.0;
                    for (size_t d = 0; d < m_mediating_pop[i].size(); ++d) {
                        double diff = m_mediating_pop[i][d] - opt.sol[d];
                        d2 += diff * diff;
                    }
                    min_dist = std::min(min_dist, std::sqrt(d2));
                }

                if (min_dist >= clear_r) {
                    if (strictLocalOptimalityCheck(m_mediating_pop[i], m_mediating_fitness[i])) {
                        // Added new: Add to tabu list when newly found
                        addTabuPeak(m_mediating_pop[i]);
                        m_found_optima.push_back({ m_mediating_pop[i], m_mediating_fitness[i], -1 });
                    }
                }
            }
            pruneFoundOptimaArchive(8);
        }


        // PER-BASIN SEPARATE REFINEMENT
        void refineFoundOptimaPerBasin() {
            clearCaches();
            if (m_found_optima.empty()) return;

            std::mt19937 rng(99999u);
            std::normal_distribution<double> noise(-0.03, 0.03);
            const int total_dim = m_bench->getNumVariables();

            for (size_t basin_idx = 0; basin_idx < m_found_optima.size(); ++basin_idx) {
                auto& opt = m_found_optima[basin_idx];

                const int mini_n = 15;
                std::vector<std::vector<double>> mini_pop(mini_n, opt.sol);
                std::vector<double> mini_fit(mini_n, opt.fitness);

                for (int i = 1; i < mini_n; ++i) {
                    for (auto& v : mini_pop[i]) {
                        v += noise(rng);
                        v = std::max(0.0, std::min(1.0, v));
                    }
                    mini_fit[i] = evalJointCached(mini_pop[i]);
                }

                for (int g = 0; g < 40; ++g) {
                    std::uniform_int_distribution<int> pick(0, mini_n - 1);
                    std::uniform_real_distribution<double> unif(0.0, 1.0);

                    for (int i = 0; i < mini_n; ++i) {
                        int r1 = i, r2 = i, r3 = i, att = 0;
                        while (r1 == i && ++att < 100) r1 = pick(rng);
                        att = 0; while ((r2 == i || r2 == r1) && ++att < 100) r2 = pick(rng);
                        att = 0; while ((r3 == i || r3 == r1 || r3 == r2) && ++att < 100) r3 = pick(rng);
                        if (r1 == i || r2 == i || r3 == i) continue;

                        std::vector<double> trial(mini_pop[i]);
                        const int jrand = static_cast<int>(unif(rng) * total_dim);
                        const double ref_F = 0.3;
                        const double ref_CR = 0.95;
                        for (int d = 0; d < total_dim; ++d) {
                            if (unif(rng) < ref_CR || d == jrand) {
                                double v = mini_pop[r1][d] + ref_F * (mini_pop[r2][d] - mini_pop[r3][d]);
                                v = std::max(0.0, std::min(1.0, v));
                                trial[d] = v;
                            }
                        }
                        double tf = evalJointCached(trial);
                        if (tf >= mini_fit[i]) {
                            mini_pop[i] = std::move(trial);
                            mini_fit[i] = tf;
                        }
                    }
                }

                auto best_it = std::max_element(mini_fit.begin(), mini_fit.end());
                size_t bi = static_cast<size_t>(std::distance(mini_fit.begin(), best_it));
                opt.sol = mini_pop[bi];
                opt.fitness = mini_fit[bi];
                opt.basin_id = static_cast<int>(basin_idx);
            }
        }


        // REFINE CANDIDATE LIST - Global refinement
        void refineCandidateList(std::vector<std::vector<double>>& candidates) {
            clearCaches();
            if (candidates.empty()) return;
            const int total_dim = m_bench->getNumVariables();
            const double REFINE_F = 0.45;
            const double REFINE_CR = 0.95;
            const int mini_n = 30;
            const int mini_g = 60;

            std::vector<double> fits(candidates.size());
            for (size_t i = 0; i < candidates.size(); ++i)
                fits[i] = evalJointCached(candidates[i]);

            std::vector<size_t> order(candidates.size());
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(),
                [&](size_t a, size_t b) { return fits[a] > fits[b]; });

            std::vector<std::vector<double>> selected;
            std::vector<double> selected_fit;
            for (size_t idx : order) {
                if (selected.size() >= 80) break;
                bool too_close = false;
                for (const auto& ex : selected) {
                    double d2 = 0.0;
                    for (size_t d = 0; d < candidates[idx].size(); ++d) {
                        double diff = candidates[idx][d] - ex[d];
                        d2 += diff * diff;
                    }
                    if (std::sqrt(d2) < 0.03) { too_close = true; break; }
                }
                if (!too_close) {
                    selected.push_back(candidates[idx]);
                    selected_fit.push_back(fits[idx]);
                }
            }
            if (selected.empty()) return;

            std::mt19937 rng(987654u);
            std::normal_distribution<double> noise(0.0, 0.015);
            std::uniform_real_distribution<double> unif(0.0, 1.0);

            for (size_t s = 0; s < selected.size(); ++s) {
                std::vector<std::vector<double>> mp(mini_n, selected[s]);
                std::vector<double> mf(mini_n, selected_fit[s]);
                for (int i = 1; i < mini_n; ++i) {
                    for (auto& v : mp[i]) {
                        v += noise(rng);
                        v = std::max(0.0, std::min(1.0, v));
                    }
                    mf[i] = evalJointCached(mp[i]);
                }

                for (int g = 0; g < mini_g; ++g) {
                    std::uniform_int_distribution<int> pick(0, mini_n - 1);
                    for (int i = 0; i < mini_n; ++i) {
                        int r1 = i, r2 = i, r3 = i, att = 0;
                        while (r1 == i && ++att < 100) r1 = pick(rng);
                        att = 0; while ((r2 == i || r2 == r1) && ++att < 100) r2 = pick(rng);
                        att = 0; while ((r3 == i || r3 == r1 || r3 == r2) && ++att < 100) r3 = pick(rng);
                        if (r1 == i || r2 == i || r3 == i) continue;

                        std::vector<double> trial(mp[i]);
                        const int jrand = static_cast<int>(unif(rng) * total_dim);
                        for (int d = 0; d < total_dim; ++d) {
                            if (unif(rng) < REFINE_CR || d == jrand) {
                                double v = mp[r1][d] + REFINE_F * (mp[r2][d] - mp[r3][d]);
                                v = std::max(0.0, std::min(1.0, v));
                                trial[d] = v;
                            }
                        }
                        double tf = evalJointCached(trial);
                        if (tf >= mf[i]) {
                            mp[i] = std::move(trial);
                            mf[i] = tf;
                        }
                    }
                }
                auto best_it = std::max_element(mf.begin(), mf.end());
                size_t bi = static_cast<size_t>(std::distance(mf.begin(), best_it));
                selected[s] = mp[bi];
                selected_fit[s] = mf[bi];
            }
            candidates = std::move(selected);
        }


        // GETTERS AND UTILITIES
        double getBestConsensusFitness() const {
            double best = -std::numeric_limits<double>::max();
            for (double f : m_mediating_fitness)
                if (std::isfinite(f) && f > best) best = f;
            return best;
        }

        std::vector<double> getBestConsensusSolution() const {
            if (m_mediating_pop.empty()) return {};
            double best = -std::numeric_limits<double>::max();
            size_t bi = 0;
            for (size_t i = 0; i < m_mediating_fitness.size(); ++i) {
                if (std::isfinite(m_mediating_fitness[i]) && m_mediating_fitness[i] > best) {
                    best = m_mediating_fitness[i];
                    bi = i;
                }
            }
            return m_mediating_pop[bi];
        }

        std::vector<std::vector<double>> getFoundOptimaSolutions() const {
            std::vector<std::vector<double>> result;
            for (const auto& opt : m_found_optima) result.push_back(opt.sol);
            return result;
        }

        std::vector<std::vector<double>> getCartesianJointCandidates() const {
            std::vector<std::vector<double>> candidates;
            const int num_dms = m_bench->getNumDMs();
            const int total_dim = m_bench->getNumVariables();
            if (num_dms == 0) return candidates;

            std::vector<int> cur_subpop(num_dms, 0);
            std::function<void(int)> dfs = [&](int dm) {
                if (dm == num_dms) {
                    std::vector<double> full_sol(total_dim, 0.5);
                    for (int d = 0; d < num_dms; ++d) {
                        const auto& idx = m_bench->getPartyIndices(d);
                        const auto& pop = m_competing_pops[d][cur_subpop[d]];
                        if (pop.empty()) return;
                        const auto& best = pop[0];
                        for (size_t k = 0; k < idx.size() && k < best.size(); ++k)
                            full_sol[idx[k]] = best[k];
                    }
                    candidates.push_back(full_sol);
                    return;
                }
                for (size_t sp = 0; sp < m_competing_pops[dm].size(); ++sp) {
                    if (m_competing_pops[dm][sp].empty()) continue;
                    cur_subpop[dm] = static_cast<int>(sp);
                    dfs(dm + 1);
                }
                };
            dfs(0);
            return candidates;
        }

        std::vector<std::vector<double>> getAllCandidateSolutions() const {
            std::vector<std::vector<double>> candidates = m_mediating_pop;
            auto refined = getFoundOptimaSolutions();
            candidates.insert(candidates.end(), refined.begin(), refined.end());
            auto cartesian = getCartesianJointCandidates();
            candidates.insert(candidates.end(), cartesian.begin(), cartesian.end());
            return candidates;
        }

        const std::vector<std::vector<double>>& getMediatingPopulation() const {
            return m_mediating_pop;
        }

        // OUTPUT FUNCTIONS
        void saveSolutionsAtGeneration(int gen_num) const {
            const std::string fn =
                m_output_directory + "/gen_" + std::to_string(gen_num) + "_solutions.txt";
            std::ofstream out(fn);
            if (!out.is_open()) return;
            out << std::fixed << std::setprecision(6)
                << "# gen=" << gen_num << "  format: x0 x1 ... xd joint_fitness\n";
            for (size_t i = 0; i < m_mediating_pop.size(); i++) {
                for (const double v : m_mediating_pop[i]) out << v << " ";
                out << m_mediating_fitness[i] << "\n";
            }
        }

        void saveConvergenceCurve() const {
            const std::string fn = m_output_directory + "/../convergence.txt";
            std::ofstream out(fn);
            if (!out.is_open()) return;
            out << "# generation best_joint_fitness\n";
            for (const auto& p : m_convergence_log)
                out << p.first << " " << std::fixed << std::setprecision(6) << p.second << "\n";
        }

        void outputPartyLandscapes(const std::string& output_dir) const {
            std::filesystem::create_directories(output_dir);
            const int res = 200;
            const int num_dms = m_bench->getNumDMs();

            std::vector<std::vector<std::vector<double>>> contexts;
            std::vector<std::string> ctx_names;

            std::vector<std::vector<double>> neutral(num_dms);
            for (int dm = 0; dm < num_dms; ++dm)
                neutral[dm].assign(m_bench->getPartyDim(dm), 0.5);
            contexts.push_back(neutral);
            ctx_names.push_back("neutral");

            if (!m_mediating_pop.empty()) {
                auto best_it = std::max_element(m_mediating_fitness.begin(), m_mediating_fitness.end());
                size_t bi = static_cast<size_t>(std::distance(m_mediating_fitness.begin(), best_it));
                std::vector<std::vector<double>> best_ctx(num_dms);
                for (int dm = 0; dm < num_dms; ++dm) {
                    const auto& idx = m_bench->getPartyIndices(dm);
                    best_ctx[dm].resize(idx.size());
                    for (size_t k = 0; k < idx.size(); ++k)
                        best_ctx[dm][k] = m_mediating_pop[bi][idx[k]];
                }
                contexts.push_back(best_ctx);
                ctx_names.push_back("global_best");
            }

            if (!m_context_archive.empty()) {
                std::vector<std::vector<double>> rand_ctx(num_dms);
                for (int dm = 0; dm < num_dms; ++dm) {
                    if (!m_context_archive[dm].empty())
                        rand_ctx[dm] = m_context_archive[dm].back();
                    else
                        rand_ctx[dm].assign(m_bench->getPartyDim(dm), 0.5);
                }
                contexts.push_back(rand_ctx);
                ctx_names.push_back("diverse");
            }

            for (int dm = 0; dm < num_dms; ++dm) {
                const auto& idx = m_bench->getPartyIndices(dm);
                if (idx.size() != 2) continue;

                int g1 = idx[0], g2 = idx[1];

                for (size_t ctx_id = 0; ctx_id < contexts.size(); ++ctx_id) {
                    std::ostringstream fname;
                    fname << output_dir << "/party" << dm
                        << "_" << ctx_names[ctx_id]
                        << "_d" << g1 << "_d" << g2 << ".txt";
                    std::ofstream f(fname.str());
                    f << std::fixed << std::setprecision(6);
                    f << "# Party " << dm << " vars (" << g1 << "," << g2
                        << ") | context=" << ctx_names[ctx_id] << "\n";

                    for (int i = 0; i < res; ++i) {
                        for (int j = 0; j < res; ++j) {
                            std::vector<double> full_sol(m_bench->getNumVariables(), 0.5);
                            full_sol[g1] = double(i) / (res - 1);
                            full_sol[g2] = double(j) / (res - 1);

                            for (int other = 0; other < num_dms; ++other) {
                                if (other == dm) continue;
                                const auto& oth_idx = m_bench->getPartyIndices(other);
                                const auto& ctx = contexts[ctx_id][other];
                                for (size_t k = 0; k < oth_idx.size(); ++k)
                                    full_sol[oth_idx[k]] = ctx[k];
                            }

                            double fit = evalJointCached(full_sol);
                            f << full_sol[g1] << " " << full_sol[g2] << " " << fit << "\n";
                        }
                    }
                }
            }
        }

    private:
        // ELITE ARCHIVE
        void buildEliteArchive(int dm,
            const std::vector<std::vector<double>>& pool,
            const std::vector<double>& pool_fit) {
            const int dm_dim = m_bench->getPartyDim(dm);

            m_elite_archive[dm].clear();
            m_elite_fitness[dm].clear();
            const double min_sep = 0.05;

            std::vector<int> order(static_cast<size_t>(pool.size()));
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(),
                [&](int a, int b) { return pool_fit[a] > pool_fit[b]; });

            for (int idx : order) {
                if (static_cast<int>(m_elite_archive[dm].size()) >= m_K) break;
                bool too_close = false;
                for (const auto& ex : m_elite_archive[dm]) {
                    double d2 = 0.0;
                    for (int d = 0; d < dm_dim; d++) {
                        double diff = pool[idx][d] - ex[d];
                        d2 += diff * diff;
                    }
                    if (std::sqrt(d2) < min_sep) { too_close = true; break; }
                }
                if (!too_close) {
                    m_elite_archive[dm].push_back(pool[idx]);
                    m_elite_fitness[dm].push_back(pool_fit[idx]);
                }
            }
        }

        void updateEliteArchive(int dm) {
            for (size_t sp = 0; sp < m_competing_pops[dm].size(); sp++) {
                if (m_competing_pops[dm][sp].empty()) continue;
                const auto& best_ind = m_competing_pops[dm][sp][0];
                const double best_f = m_competing_fitness[dm][sp][0];

                if (m_elite_archive[dm].size() < static_cast<size_t>(m_K)) {
                    m_elite_archive[dm].push_back(best_ind);
                    m_elite_fitness[dm].push_back(best_f);
                }
                else {
                    auto worst_it = std::min_element(m_elite_fitness[dm].begin(),
                        m_elite_fitness[dm].end());
                    if (best_f > *worst_it) {
                        const size_t wi = static_cast<size_t>(
                            std::distance(m_elite_fitness[dm].begin(), worst_it));
                        m_elite_archive[dm][wi] = best_ind;
                        m_elite_fitness[dm][wi] = best_f;
                    }
                }
            }
        }

        // CONSENSUS CANDIDATE INJECTION
        void injectConsensusCandidates() {
            const int total_dim = m_bench->getNumVariables();
            const int num_dms = m_bench->getNumDMs();

            // Best context per party
            {
                std::vector<double> cand(total_dim, 0.5);
                for (int dm = 0; dm < num_dms; dm++) {
                    const auto& idx = m_bench->getPartyIndices(dm);
                    if (!m_context_archive[dm].empty()) {
                        const auto& ctx = m_context_archive[dm][0];
                        for (size_t i = 0; i < idx.size() && i < ctx.size(); i++)
                            cand[idx[i]] = ctx[i];
                    }
                }
                insertIfBetter(cand);
            }

            // Cartesian product of TOP-3 contexts per party
            const int TOP = 3;
            std::vector<std::vector<int>> ranges(num_dms);
            for (int dm = 0; dm < num_dms; ++dm) {
                int n = std::min(TOP, static_cast<int>(m_context_archive[dm].size()));
                ranges[dm].resize(n);
                std::iota(ranges[dm].begin(), ranges[dm].end(), 0);
            }

            std::vector<int> cur(num_dms, 0);
            std::function<void(int)> dfs = [&](int dm) {
                if (dm == num_dms) {
                    std::vector<double> cand(total_dim, 0.5);
                    for (int d = 0; d < num_dms; ++d) {
                        const auto& idx = m_bench->getPartyIndices(d);
                        if (cur[d] < static_cast<int>(m_context_archive[d].size())) {
                            const auto& ctx = m_context_archive[d][cur[d]];
                            for (size_t i = 0; i < idx.size() && i < ctx.size(); i++)
                                cand[idx[i]] = ctx[i];
                        }
                    }
                    insertIfBetter(cand);
                    return;
                }
                for (int i : ranges[dm]) {
                    cur[dm] = i;
                    dfs(dm + 1);
                }
                };
            dfs(0);

            // Random elite mixes
            std::mt19937 rng(static_cast<unsigned>(m_generation * 77777u));
            for (int combo = 0; combo < 10; ++combo) {
                std::vector<double> cand(total_dim, 0.5);
                for (int dm = 0; dm < num_dms; ++dm) {
                    if (m_context_archive[dm].empty()) continue;
                    const auto& idx = m_bench->getPartyIndices(dm);
                    std::uniform_int_distribution<int> pick(0, static_cast<int>(m_context_archive[dm].size()) - 1);
                    const auto& ctx = m_context_archive[dm][pick(rng)];
                    for (size_t i = 0; i < idx.size() && i < ctx.size(); i++)
                        cand[idx[i]] = ctx[i];
                }
                insertIfBetter(cand);
            }

            // Random subspace combinations
            std::mt19937 rng_wild(static_cast<unsigned>(m_generation * 88888u));
            for (int wild = 0; wild < 5; ++wild) {
                std::vector<double> cand(total_dim, 0.5);
                for (int dm = 0; dm < num_dms; ++dm) {
                    const auto& idx = m_bench->getPartyIndices(dm);
                    if (m_competing_pops[dm].empty() || m_context_archive[dm].empty()) continue;

                    std::uniform_int_distribution<int> pick_sp(0, static_cast<int>(m_competing_pops[dm].size()) - 1);
                    std::uniform_int_distribution<int> pick_ctx(0, static_cast<int>(m_context_archive[dm].size()) - 1);
                    std::uniform_real_distribution<double> coin(0.0, 1.0);

                    const auto& src = (coin(rng_wild) < 0.5 && !m_competing_pops[dm][pick_sp(rng_wild)].empty())
                        ? m_competing_pops[dm][pick_sp(rng_wild)][0]
                        : m_context_archive[dm][pick_ctx(rng_wild)];

                    for (size_t i = 0; i < idx.size() && i < src.size(); ++i)
                        cand[idx[i]] = src[i];
                }
                insertIfBetter(cand);
            }
        }

        void injectCartesianFromSubPops() {
            const int num_dms = m_bench->getNumDMs();
            const int total_dim = m_bench->getNumVariables();
            std::vector<int> cur_subpop(num_dms, 0);
            std::function<void(int)> dfs = [&](int dm) {
                if (dm == num_dms) {
                    std::vector<double> full_sol(total_dim, 0.5);
                    for (int d = 0; d < num_dms; ++d) {
                        const auto& idx = m_bench->getPartyIndices(d);
                        const auto& pop = m_competing_pops[d][cur_subpop[d]];
                        if (pop.empty()) return;
                        const auto& best = pop[0];
                        for (size_t k = 0; k < idx.size() && k < best.size(); ++k)
                            full_sol[idx[k]] = best[k];
                    }
                    insertIfBetter(full_sol);
                    return;
                }
                for (size_t sp = 0; sp < m_competing_pops[dm].size(); ++sp) {
                    if (m_competing_pops[dm][sp].empty()) continue;
                    cur_subpop[dm] = static_cast<int>(sp);
                    dfs(dm + 1);
                }
                };
            dfs(0);
        }

        void insertIfBetter(const std::vector<double>& cand) {
            const double raw_cf = evalJointCached(cand);
            if (!std::isfinite(raw_cf)) return;

            // Added new: Reject tabu candidates
            if (isTabu(cand)) return;

            for (size_t i = 0; i < m_mediating_pop.size(); ++i) {
                if (!std::isfinite(m_mediating_fitness[i])) continue;
                double d2 = 0.0;
                for (size_t d = 0; d < cand.size(); ++d) {
                    double diff = cand[d] - m_mediating_pop[i][d];
                    d2 += diff * diff;
                }
                if (std::sqrt(d2) < 0.03 && raw_cf < m_mediating_fitness[i] + 0.5) {
                    return;
                }
            }

            double worst_f = std::numeric_limits<double>::max();
            size_t wi = 0;
            for (size_t i = 0; i < m_mediating_fitness.size(); ++i) {
                double f = std::isfinite(m_mediating_fitness[i])
                    ? m_mediating_fitness[i]
                    : -std::numeric_limits<double>::max();
                if (f < worst_f) { worst_f = f; wi = i; }
            }
            if (raw_cf > worst_f) {
                m_mediating_pop[wi] = cand;
                m_mediating_fitness[wi] = raw_cf;
            }
        }

        void forceInsertMediating(const std::vector<double>& sol, double fit) {
            double worst_f = std::numeric_limits<double>::max();
            size_t wi = 0;
            for (size_t i = 0; i < m_mediating_fitness.size(); ++i) {
                double f = std::isfinite(m_mediating_fitness[i]) ? m_mediating_fitness[i] : -1e30;
                if (f < worst_f) { worst_f = f; wi = i; }
            }
            m_mediating_pop[wi] = sol;
            m_mediating_fitness[wi] = fit;
        }

        // CONTEXT ARCHIVE PRUNING
        void pruneContextArchiveDiverse(int dm, double min_sep) {
            auto& arch = m_context_archive[dm];
            auto& fit = m_context_fitness[dm];
            if (arch.empty()) return;

            std::vector<size_t> order(arch.size());
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(),
                [&](size_t a, size_t b) { return fit[a] > fit[b]; });

            std::vector<std::vector<double>> new_arch;
            std::vector<double> new_fit;

            for (size_t idx : order) {
                if (new_arch.size() >= static_cast<size_t>(m_K)) break;
                bool too_close = false;
                for (const auto& ex : new_arch) {
                    double d2 = 0.0;
                    for (size_t d = 0; d < arch[idx].size(); ++d) {
                        double diff = arch[idx][d] - ex[d];
                        d2 += diff * diff;
                    }
                    if (std::sqrt(d2) < min_sep) { too_close = true; break; }
                }
                if (!too_close) {
                    new_arch.push_back(arch[idx]);
                    new_fit.push_back(fit[idx]);
                }
            }
            arch = std::move(new_arch);
            fit = std::move(new_fit);
        }

        // FOUND OPTIMA PRUNING
        void pruneFoundOptimaArchive(size_t hard_cap = 50) {
            if (m_found_optima.size() <= 1) return;

            std::vector<size_t> order(m_found_optima.size());
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(),
                [&](size_t a, size_t b) { return m_found_optima[a].fitness > m_found_optima[b].fitness; });

            std::vector < FoundOptimum > pruned;
            for (size_t idx : order) {
                bool too_close = false;
                // Added new: Use adaptive clearing radius
                double clear_r = getAdaptiveClearingRadius();

                for (const auto& ex : pruned) {
                    double d2 = 0.0;
                    for (size_t d = 0; d < m_found_optima[idx].sol.size(); ++d) {
                        double diff = m_found_optima[idx].sol[d] - ex.sol[d];
                        d2 += diff * diff;
                    }
                    if (std::sqrt(d2) < clear_r) { too_close = true; break; }
                }
                if (!too_close) {
                    pruned.push_back(m_found_optima[idx]);
                    if (pruned.size() >= hard_cap) break;
                }
            }
            m_found_optima = std::move(pruned);
        }

        // MEDIATING POPULATION MAINTENANCE
        void clearMediatingPopulation() {
            std::mt19937 rng(static_cast<unsigned>(m_generation * 1234567u));
            std::uniform_real_distribution<double> unif(0.0, 1.0);

            const int n_clear = std::max(1, static_cast<int>(m_mediating_pop.size()) / 10);
            for (int c = 0; c < n_clear; ++c) {
                size_t wi = 0;
                double min_nn_dist = std::numeric_limits<double>::max();

                double best_f = -1e30;
                size_t best_i = 0;
                for (size_t i = 0; i < m_mediating_fitness.size(); ++i) {
                    if (std::isfinite(m_mediating_fitness[i]) && m_mediating_fitness[i] > best_f) {
                        best_f = m_mediating_fitness[i]; best_i = i;
                    }
                }

                for (size_t i = 0; i < m_mediating_pop.size(); ++i) {
                    if (i == best_i) continue;
                    double nearest = std::numeric_limits<double>::max();
                    for (size_t j = 0; j < m_mediating_pop.size(); ++j) {
                        if (i == j) continue;
                        double d2 = 0.0;
                        for (size_t d = 0; d < m_mediating_pop[i].size(); ++d) {
                            double diff = m_mediating_pop[i][d] - m_mediating_pop[j][d];
                            d2 += diff * diff;
                        }
                        nearest = std::min(nearest, std::sqrt(d2));
                    }
                    if (nearest < min_nn_dist) {
                        min_nn_dist = nearest;
                        wi = i;
                    }
                }

                std::vector<double> rsol(m_mediating_pop[wi].size());
                bool valid = false;
                int attempts = 0;
                while (!valid && attempts < 100) {
                    for (auto& v : rsol) v = unif(rng);
                    valid = true;
                    // Added new :Tabu check for replacement individuals
                    if (isTabu(rsol)) valid = false;
                    for (const auto& opt : m_found_optima) {
                        if (opt.fitness < 60.0) continue;
                        double d2 = 0.0;
                        for (size_t d = 0; d < rsol.size(); ++d) {
                            double diff = rsol[d] - opt.sol[d];
                            d2 += diff * diff;
                        }
                        if (std::sqrt(d2) < 0.10) { valid = false; break; }
                    }
                    ++attempts;
                }
                m_mediating_pop[wi] = rsol;
                m_mediating_fitness[wi] = evalJointCached(rsol);
            }
        }

        void mergeContextsFromMediating() {
            const int num_dms = m_bench->getNumDMs();
            for (int dm = 0; dm < num_dms; ++dm) {
                const auto& idx = m_bench->getPartyIndices(dm);
                const int dm_dim = m_bench->getPartyDim(dm);

                std::vector<size_t> order(m_mediating_pop.size());
                std::iota(order.begin(), order.end(), 0);
                std::sort(order.begin(), order.end(),
                    [&](size_t a, size_t b) { return m_mediating_fitness[a] > m_mediating_fitness[b]; });

                for (int r = 0; r < 2 && r < static_cast<int>(order.size()); ++r) {
                    std::vector<double> ctx(dm_dim);
                    for (int d = 0; d < dm_dim; ++d) ctx[d] = m_mediating_pop[order[r]][idx[d]];
                    m_context_archive[dm].push_back(std::move(ctx));
                    m_context_fitness[dm].push_back(m_mediating_fitness[order[r]]);
                }

                std::mt19937 rng_div(static_cast<unsigned>(m_generation * 44444u + dm));
                std::uniform_int_distribution<int> pick_rand(0, static_cast<int>(m_mediating_pop.size()) - 1);
                for (int rand_i = 0; rand_i < 5; ++rand_i) {
                    int ri = pick_rand(rng_div);
                    std::vector<double> ctx(dm_dim);
                    for (int d = 0; d < dm_dim; ++d) ctx[d] = m_mediating_pop[ri][idx[d]];
                    m_context_archive[dm].push_back(std::move(ctx));
                    m_context_fitness[dm].push_back(m_mediating_fitness[ri]);
                }

                {
                    std::vector<int> candidates;
                    for (int i = 0; i < static_cast<int>(m_mediating_pop.size()); ++i) candidates.push_back(i);
                    std::shuffle(candidates.begin(), candidates.end(), rng_div);
                    const int n_diverse = 5;
                    for (int d = 0; d < n_diverse && !candidates.empty(); ++d) {
                        double best_min_dist = -1.0;
                        int best_idx = -1;
                        for (int ci : candidates) {
                            double min_dist = 1e30;
                            for (const auto& arch : m_context_archive[dm]) {
                                double d2 = 0.0;
                                for (int dd = 0; dd < dm_dim; ++dd) {
                                    double diff = m_mediating_pop[ci][idx[dd]] - arch[dd];
                                    d2 += diff * diff;
                                }
                                min_dist = std::min(min_dist, std::sqrt(d2));
                            }
                            if (min_dist > best_min_dist) {
                                best_min_dist = min_dist;
                                best_idx = ci;
                            }
                        }
                        if (best_idx >= 0) {
                            std::vector<double> ctx(dm_dim);
                            for (int dd = 0; dd < dm_dim; ++dd) ctx[dd] = m_mediating_pop[best_idx][idx[dd]];
                            m_context_archive[dm].push_back(std::move(ctx));
                            m_context_fitness[dm].push_back(m_mediating_fitness[best_idx]);
                            candidates.erase(std::remove(candidates.begin(), candidates.end(), best_idx), candidates.end());
                        }
                    }
                }

                if (!m_found_optima.empty()) {
                    for (const auto& opt : m_found_optima) {
                        double best_d2 = 1e30;
                        size_t best_i = 0;
                        for (size_t i = 0; i < m_mediating_pop.size(); ++i) {
                            double d2 = 0.0;
                            for (size_t d = 0; d < m_mediating_pop[i].size(); ++d) {
                                double diff = m_mediating_pop[i][d] - opt.sol[d];
                                d2 += diff * diff;
                            }
                            if (d2 < best_d2) { best_d2 = d2; best_i = i; }
                        }
                        std::vector<double> ctx(dm_dim);
                        for (int d = 0; d < dm_dim; ++d) ctx[d] = m_mediating_pop[best_i][idx[d]];
                        m_context_archive[dm].push_back(std::move(ctx));
                        m_context_fitness[dm].push_back(m_mediating_fitness[best_i]);
                    }
                }

                std::uniform_real_distribution<double> pure_unif(0.0, 1.0);
                for (int pure = 0; pure < 5; ++pure) {
                    std::vector<double> ctx(dm_dim);
                    for (int d = 0; d < dm_dim; ++d) ctx[d] = pure_unif(rng_div);
                    m_context_archive[dm].push_back(std::move(ctx));
                    m_context_fitness[dm].push_back(-1e30);
                }

                pruneContextArchiveDiverse(dm, 0.10);
            }
        }
    };


    // CEC-2015 METRICS
    struct CEC2015Metrics {
        double SR = 0.0;
        double ANOF = 0.0;
        double SP = 0.0;
        double MPR = 0.0;
    };

    CEC2015Metrics evaluateCEC2015Metrics(
        const std::vector<std::vector<double>>& candidates,
        FreePeaks* problem,
        Environment* env,
        double eps_x = 0.05,
        double eps_f = 2.0)
    {
        CEC2015Metrics m;
        const auto* opt = problem->optima();
        if (!opt || opt->numberSolutions() == 0) return m;

        const int n_opt = static_cast<int>(opt->numberSolutions());
        std::vector<bool> found(static_cast<size_t>(n_opt), false);
        int total_found = 0;

        for (const auto& cand : candidates) {
            std::unique_ptr<SolutionBase> sol(problem->createSolution());
            auto& sv = dynamic_cast<Solution < VariableVector<Real>>& > (*sol)
                .variable().vector();
            for (size_t d = 0; d < cand.size() && d < sv.size(); d++)
                sv[d] = static_cast<Real>(cand[d]);
            sol->evaluate(env, false);
            const double cf = sol->objective(0);
            if (!std::isfinite(cf)) continue;

            for (int i = 0; i < n_opt; i++) {
                if (found[static_cast<size_t>(i)]) continue;
                const auto& ov = opt->solution(i).variable().vector();
                double d2 = 0.0;
                for (size_t d = 0; d < cand.size() && d < ov.size(); d++) {
                    double diff = cand[d] - static_cast<double>(ov[d]);
                    d2 += diff * diff;
                }
                const double of = opt->solution(i).objective(0);
                if (std::sqrt(d2) < eps_x && std::abs(cf - of) < eps_f) {
                    found[static_cast<size_t>(i)] = true;
                    ++total_found;
                }
            }
        }

        m.ANOF = static_cast<double>(total_found) / static_cast<double>(n_opt);
        m.SR = (total_found == n_opt) ? 100.0 : 0.0;
        m.SP = (m.SR > 0.0)
            ? static_cast<double>(candidates.size()) * n_opt
            : std::numeric_limits<double>::infinity();

        double global_f = -std::numeric_limits<double>::max();
        for (int i = 0; i < n_opt; i++)
            global_f = std::max(global_f, opt->solution(i).objective(0));

        double sum_num = 0.0, sum_den = 0.0;
        for (int i = 0; i < n_opt; i++) {
            const double w = (global_f - opt->solution(i).objective(0)) + 1.0;
            if (found[static_cast<size_t>(i)]) sum_num += w;
            sum_den += w;
        }
        m.MPR = (sum_den > 0.0) ? (sum_num / sum_den) : 0.0;
        return m;
    }


    // MANUAL PEAK CENTER INSPECTION
    void inspectPeakCenters(FreePeaks* problem, Environment* env) {
        const auto* opt = problem->optima();
        if (!opt || opt->numberSolutions() == 0) {
            std::cout << "[INSPECT] No optima recorded in problem.\n";
            return;
        }

        std::cout << "\n>>> MANUAL PEAK CENTER INSPECTION\n";
        std::cout << "    Total peaks: " << opt->numberSolutions() << "\n";
        std::cout << "    Format: [x0, x1, x2, x3] -> fitness\n";

        double max_fitness = -std::numeric_limits<double>::max();
        int best_idx = -1;

        for (int i = 0; i < static_cast<int>(opt->numberSolutions()); ++i) {
            const auto& sol = opt->solution(i);
            const auto& vars = sol.variable().vector();
            const double f = sol.objective(0);

            std::cout << "    Peak " << (i + 1) << ": [";
            for (size_t d = 0; d < vars.size(); ++d) {
                std::cout << std::fixed << std::setprecision(4) << static_cast<double>(vars[d]);
                if (d + 1 < vars.size()) std::cout << ", ";
            }
            std::cout << "] -> " << std::fixed << std::setprecision(4) << f << "\n";

            if (f > max_fitness) {
                max_fitness = f;
                best_idx = i;
            }
        }

        std::cout << "    Global optimum: Peak " << (best_idx + 1)
            << " with fitness " << std::fixed << std::setprecision(4) << max_fitness << "\n";

        std::cout << "\n    Direct evaluation verification:\n";
        for (int i = 0; i < static_cast<int>(opt->numberSolutions()); ++i) {
            const auto& sol = opt->solution(i);
            const auto& vars = sol.variable().vector();

            std::vector<double> full_sol(vars.size());
            for (size_t d = 0; d < vars.size(); ++d)
                full_sol[d] = static_cast<double>(vars[d]);

            std::unique_ptr<SolutionBase> test_sol(problem->createSolution());
            auto& sv = dynamic_cast<Solution < VariableVector<Real>>& > (*test_sol)
                .variable().vector();
            for (size_t d = 0; d < full_sol.size() && d < sv.size(); ++d)
                sv[d] = static_cast<Real>(full_sol[d]);
            test_sol->evaluate(env, false);
            const double eval_f = test_sol->objective(0);

            std::cout << "    Peak " << (i + 1) << " direct eval: "
                << std::fixed << std::setprecision(4) << eval_f;
            if (std::abs(eval_f - sol.objective(0)) > 0.01) {
                std::cout << " [MISMATCH with stored: " << sol.objective(0) << "]";
            }
            std::cout << "\n";
        }
    }


    // PROBLEM CREATION
    void createComplex4D2PartyProblem(const std::string& problem_name = "complex_4d_2p") {
        using namespace ofec;
        using namespace free_peaks;

        const std::string problem_dir = "multiparty_multimodal/";
        const std::string sub_dir = "sop/";

        std::filesystem::create_directories(
            g_working_directory + "instance/problem/continuous/free_peaks/" + problem_dir);
        std::filesystem::create_directories(
            g_working_directory + "instance/problem/continuous/free_peaks/" + sub_dir);
        std::filesystem::create_directories(
            g_working_directory + "instance/problem/continuous/free_peaks/subproblem/" + sub_dir);
        std::filesystem::create_directories(
            g_working_directory + "instance/problem/continuous/free_peaks/subproblem/function/one_peak/" + sub_dir);

        const int numDim = 4;
        const int numObj = 1;
        const int numCon = 0;

        std::shared_ptr < ofec::Random > rnd(new Random(0.123));

        FreePeaks::registerFP();
        const std::string fp_name = "free_peaks";
        std::shared_ptr < Environment > env(generateEnvironmentByFactory(fp_name));
        env->recordInputParameters();
        env->initialize();
        env->setProblem(generateProblemByFactory(fp_name));
        auto* freepeak = dynamic_cast < FreePeaks* > (env->problem());

        ParameterMap fp_pm;
        fp_pm["generation_type"] = std::string("assigned");
        fp_pm["dataFile1"] = problem_dir + problem_name + ".txt";
        freepeak->inputParameters().input(fp_pm);
        freepeak->initialize(env.get());
        freepeak->setRandom(rnd);
        freepeak->setSizes(numDim, numObj, numCon);

        freepeak->setKDtree({
            { "root", {
                { "peak1", 0.20 },
                { "peak2", 0.16 },
                { "peak3", 0.13 },
                { "peak4", 0.12 },
                { "peak5", 0.11 },
                { "peak6", 0.10 },
                { "peak7", 0.10 },
                { "peak8", 0.08 }
            }}
            });

        struct PeakSpec {
            std::vector<double> center;
            double height;
            std::string shape;
            double condition;
            bool use_asym;
            int bias_party;
            double magnitude;
            double rotation;
            double pb_cond;
        };

        const std::vector < PeakSpec > peaks = {
            { {  0,   0,   0,   0  }, 90.0, "s1",   1.0,  false, 0,  0.0, 0.0,   1.0 },
            { { 40,  -30,  35, -25  }, 75.0, "s2",  10.0,  true,  0, 15.0, 0.4,  10.0 },
            { {-50,   40, -45,  35  }, 65.0, "s7",  10.0,  false, 1, 12.0, 0.6,  10.0 },
            { { 25,   25,  25,  25  }, 55.0, "s3",  50.0,  true,  0, 20.0, 0.3,  50.0 },
            { {-30,  -30, -30, -30  }, 50.0, "s5",   1.0,  false, 1,  8.0, 0.0,   1.0 },
            { { 60,  -50,  55, -45  }, 42.0, "s2", 100.0,  true,  1, 18.0, 0.5, 100.0 },
            { {-20,   60, -15,  55  }, 35.0, "s4",   1.0,  false, 0,  5.0, 0.0,   1.0 },
            { { 80,  80,  75,  80 }, 80.0, "s7", 100.0, true, 1, 12.0, 0.5, 100.0 },
        };

        const std::vector<std::string> sp_names = {
            "peak1","peak2","peak3","peak4","peak5","peak6","peak7","peak8"
        };

        for (size_t s = 0; s < sp_names.size(); s++) {
            const std::string& sname = sp_names[s];
            const PeakSpec& pk = peaks[s];

            ParameterMap spm;
            spm["subspace"] = sname;
            spm["generation_type"] = std::string("assigned");
            spm["dataFile1"] = sub_dir + problem_name + "_" + sname + ".txt";

            auto subpro(Subproblem::create());
            subpro->initialize(spm, freepeak);

            {
                auto dis(FactoryFP < DistanceBase > ::produce("Euclidean"));
                ParameterMap dp;
                dis->initialize(freepeak, sname, dp);
                subpro->setDistance(dis);
            }

            {
                ParameterMap fun_p;
                fun_p["generation_type"] = std::string("assigned");
                fun_p["dataFile1"] = sub_dir + problem_name + "_" + sname + "_fn.txt";

                auto func(FactoryFP<FunctionBase>::produce("one_peak"));
                func->initialize(freepeak, sname, fun_p);
                auto* opf = dynamic_cast < OnePeakFunction* > (func);

                auto onepeak(FactoryFP < OnePeakBase > ::produce(pk.shape));
                ParameterMap op;
                op["center_type"] = std::string("assigned");
                op["height"] = static_cast<Real>(pk.height);
                std::vector<Real> center_r(pk.center.size());
                for (size_t d = 0; d < pk.center.size(); d++) {
                    center_r[d] = static_cast<Real>(
                        std::max(-100.0, std::min(100.0, pk.center[d])));
                }
                op["center_postion"] = center_r;
                onepeak->initialize(freepeak, sname, op);
                opf->addOnePeaks(onepeak);

                subpro->setFunction(func);
            }

            {
                auto t(FactoryFP<X_TransformBase>::produce("MapXPartyBias"));
                ParameterMap tp;
                tp["party_id"] = static_cast<int>(pk.bias_party);
                tp["magnitude"] = static_cast<Real>(pk.magnitude);
                tp["rotation_angle"] = static_cast<Real>(pk.rotation);
                tp["condition_number"] = static_cast<Real>(pk.pb_cond);
                t->initialize(freepeak, sname, tp);
                subpro->addVariableTransform(t);
            }
            {
                auto t(FactoryFP<X_TransformBase>::produce("MapXIllConditioning"));
                ParameterMap tp;
                tp["condition"] = static_cast<Real>(pk.condition);
                t->initialize(freepeak, sname, tp);
                subpro->addVariableTransform(t);
            }
            {
                auto t(FactoryFP<X_TransformBase>::produce("MapXIrregularity"));
                ParameterMap tp;
                t->initialize(freepeak, sname, tp);
                subpro->addVariableTransform(t);
            }
            if (pk.use_asym) {
                auto t(FactoryFP<X_TransformBase>::produce("MapXAssymetrix"));
                ParameterMap tp;
                tp["alpha"] = static_cast<Real>(0.2 + 0.05 * static_cast<double>(s));
                t->initialize(freepeak, sname, tp);
                subpro->addVariableTransform(t);
            }
            {
                auto t(FactoryFP<X_TransformBase>::produce("MapXDeceptive"));
                ParameterMap tp;
                tp["alpha"] = static_cast<Real>(0.05 + 0.02 * static_cast<double>(s));
                t->initialize(freepeak, sname, tp);
                subpro->addVariableTransform(t);
            }
            {
                auto t(FactoryFP<X_TransformBase>::produce("MapXLinkage"));
                ParameterMap tp;
                tp["beta"] = static_cast<Real>(0.1);
                t->initialize(freepeak, sname, tp);
                subpro->addVariableTransform(t);
            }
            {
                auto ot(FactoryFP < Y_TransformBase > ::produce("map_objective"));
                ParameterMap otp;
                ot->initialize(freepeak, sname, otp);
                subpro->addObjectiveTransform(ot);
            }

            freepeak->setSubproblem(sname, subpro);
        }

        freepeak->bindData();

        freepeak->inputParameters().at("generation_type")->setValue("read_file");
        for (auto& it : freepeak->subspaceTree().name_box_subproblem)
            if (it.second.second) {
                it.second.second->inputParameters().at("generation_type")->setValue("read_file");
                it.second.second->function()->inputParameters().at("generation_type")->setValue("read_file");
            }
        freepeak->recordInputParameters();
        freepeak->outputTotalFile();

        std::cout << ">>> Problem '" << problem_name
            << "' created: " << peaks.size()
            << " peaks, dim=" << numDim << ", 2 parties\n";
    }


    // LANDSCAPE SLICE OUTPUT
    void outputLandscapeSlices(Environment* env, const std::string& output_dir) {
        std::filesystem::create_directories(output_dir);
        const int res = 200;

        auto* con_pro = CAST_CONOP(env->problem());
        const int dim = static_cast<int>(con_pro->numberVariables());

        std::cout << ">>> Sampling landscape slices (res=" << res << ")..." << std::flush;

        for (int d1 = 0; d1 < dim; d1++) {
            for (int d2 = d1 + 1; d2 < dim; d2++) {
                std::ostringstream fname;
                fname << output_dir << "/slice_d" << d1 << "_d" << d2 << ".txt";
                std::ofstream f(fname.str());
                f << std::fixed << std::setprecision(6);
                f << "# x" << d1 << " x" << d2 << " f(x) [all others=0.5]\n";

                for (int i = 0; i < res; i++) {
                    for (int j = 0; j < res; j++) {
                        auto sol = env->problem()->createSolution();
                        auto& v = dynamic_cast<Solution < VariableVector<Real>>& > (*sol)
                            .variable().vector();
                        v.assign(static_cast<size_t>(dim), Real(0.5));
                        v[static_cast<size_t>(d1)] = Real(i) / Real(res - 1);
                        v[static_cast<size_t>(d2)] = Real(j) / Real(res - 1);
                        sol->evaluate(env, false);
                        f << (double(i) / (res - 1)) << " "
                            << (double(j) / (res - 1)) << " "
                            << sol->objective(0) << "\n";
                    }
                }
            }
        }
        std::cout << " done.\n";
    }


    // MAIN EXPERIMENT RUNNER
    void runMPMCoEAExperiment() {
        ofec::registerInstance();

        ofec::g_working_directory =
            R"(E:\HITSZ\Research\Multimodal_Multiparty_Optimization\ThesisProject\Data\ofec_data_new\)";
        for (auto& c : ofec::g_working_directory)
            if (c == '\\') c = '/';

        const std::string vis_base =
            R"(E:\HITSZ\Research\Multimodal_Multiparty_Optimization\ThesisProject\Visualization)";

        const std::string problem_name = "complex_4d_2p";
        std::cout << "\n>>> [MPM-CoEA v2.0] Creating 4D 2-party benchmark...\n";
        createComplex4D2PartyProblem(problem_name);

        const std::string fp_name = "free_peaks";
        auto* env_raw = generateEnvironmentByFactory(fp_name);
        env_raw->recordInputParameters();
        env_raw->initialize();
        env_raw->setProblem(generateProblemByFactory(fp_name));

        ParameterMap load_pm;
        load_pm["generation_type"] = std::string("read_file");
        load_pm["dataFile1"] = std::string("multiparty_multimodal/") + problem_name + ".txt";
        env_raw->problem()->inputParameters().input(load_pm);
        env_raw->problem()->recordInputParameters();
        env_raw->initializeProblem(0.5);

        auto* fp = CAST_FPs(env_raw->problem());
        if (!fp) throw std::runtime_error("Problem is not FreePeaks.");

        const auto* opt = fp->optima();
        const int n_opt = opt ? static_cast<int>(opt->numberSolutions()) : 0;
        std::cout << ">>> Loaded: dim=" << fp->numberVariables()
            << " | parties=2 | num_optima=" << n_opt << "\n";

        inspectPeakCenters(fp, env_raw);

        outputLandscapeSlices(env_raw, vis_base + "/landscape_before");

        auto bench = std::make_shared < MPMMO_Benchmark > (fp, 2, env_raw);
        auto env_sp = std::shared_ptr < Environment > (env_raw, [](Environment*) {});

        // Tuned parameters
        MPM_CoEA solver(bench, env_sp,
            /*pop_size*/ 60,
            /*med_pop*/  600,
            /*K*/        10,
            vis_base + "/solutions");

        solver.initialize();
        solver.outputPartyLandscapes(vis_base + "/party_landscapes_before");

        solver.saveSolutionsAtGeneration(0);
        std::cout << "Gen   0 | " << std::fixed << std::setprecision(4)
            << solver.getBestConsensusFitness() << "\n";

        const int MAX_GEN = 600;

        for (int gen = 1; gen <= MAX_GEN; gen++) {
            solver.run_one_generation();

            if (gen % 5 == 0 || gen == MAX_GEN) {
                std::cout << "Gen " << std::setw(4) << gen << " | "
                    << std::fixed << std::setprecision(4)
                    << solver.getBestConsensusFitness() << "\n";
            }

            if (gen % 5 == 0 || gen == MAX_GEN)
                solver.saveSolutionsAtGeneration(gen);
        }

        solver.outputPartyLandscapes(vis_base + "/party_landscapes_after");
        solver.saveConvergenceCurve();

        solver.refineFoundOptimaPerBasin();

        auto candidates = solver.getAllCandidateSolutions();
        solver.refineCandidateList(candidates);

        // Second-pass refinement
        {
            const int total_dim = fp->numberVariables();
            const int pass2_n = 40;
            const int pass2_g = 100;
            const double P2_F = 0.30;
            const double P2_CR = 0.98;
            std::mt19937 rng(111111u);
            std::normal_distribution<double> noise(0.0, 0.01);
            std::uniform_real_distribution<double> unif(0.0, 1.0);

            std::vector<double> cfit(candidates.size());
            for (size_t i = 0; i < candidates.size(); ++i)
                cfit[i] = bench->evaluateJoint(candidates[i], env_raw);

            std::vector<size_t> order(candidates.size());
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(),
                [&](size_t a, size_t b) { return cfit[a] > cfit[b]; });
            std::vector<std::vector<double>> selected;
            for (size_t idx : order) {
                if (selected.size() >= 40) break;
                bool too_close = false;
                for (const auto& ex : selected) {
                    double d2 = 0.0;
                    for (size_t d = 0; d < candidates[idx].size(); ++d) {
                        double diff = candidates[idx][d] - ex[d];
                        d2 += diff * diff;
                    }
                    if (std::sqrt(d2) < 0.02) { too_close = true; break; }
                }
                if (!too_close) selected.push_back(candidates[idx]);
            }

            for (auto& sel : selected) {
                std::vector<std::vector<double>> mp(pass2_n, sel);
                std::vector<double> mf(pass2_n);
                for (int i = 0; i < pass2_n; ++i) {
                    if (i > 0) {
                        for (auto& v : mp[i]) {
                            v += noise(rng);
                            v = std::max(0.0, std::min(1.0, v));
                        }
                    }
                    mf[i] = bench->evaluateJoint(mp[i], env_raw);
                }
                for (int g = 0; g < pass2_g; ++g) {
                    std::uniform_int_distribution<int> pick(0, pass2_n - 1);
                    for (int i = 0; i < pass2_n; ++i) {
                        int r1 = i, r2 = i, r3 = i, att = 0;
                        while (r1 == i && ++att < 100) r1 = pick(rng);
                        att = 0; while ((r2 == i || r2 == r1) && ++att < 100) r2 = pick(rng);
                        att = 0; while ((r3 == i || r3 == r1 || r3 == r2) && ++att < 100) r3 = pick(rng);
                        if (r1 == i || r2 == i || r3 == i) continue;
                        std::vector<double> trial(mp[i]);
                        int jrand = static_cast<int>(unif(rng) * total_dim);
                        for (int d = 0; d < total_dim; ++d) {
                            if (unif(rng) < P2_CR || d == jrand) {
                                double v = mp[r1][d] + P2_F * (mp[r2][d] - mp[r3][d]);
                                v = std::max(0.0, std::min(1.0, v));
                                trial[d] = v;
                            }
                        }
                        double tf = bench->evaluateJoint(trial, env_raw);
                        if (tf >= mf[i]) { mp[i] = std::move(trial); mf[i] = tf; }
                    }
                }
                auto best_it = std::max_element(mf.begin(), mf.end());
                size_t bi = static_cast<size_t>(std::distance(mf.begin(), best_it));
                sel = mp[bi];
            }
            candidates.insert(candidates.end(), selected.begin(), selected.end());
        }

        std::cout << "\n>>> CEC-2015 Metrics (STRICT: eps_x=0.01, eps_f=1.0)\n";
        auto metrics_strict = evaluateCEC2015Metrics(candidates, fp, env_raw, 0.01, 1.0);
        std::cout << "    SR       : " << metrics_strict.SR << "%\n"
            << "    ANOF     : " << metrics_strict.ANOF << "  ("
            << static_cast<int>(std::round(metrics_strict.ANOF * n_opt))
            << "/" << n_opt << " found)\n"
            << "    MPR      : " << metrics_strict.MPR << "\n"
            << "    SP       : " << metrics_strict.SP << "\n";

        std::cout << "\n>>> CEC-2015 Metrics (LOOSE: eps_x=0.05, eps_f=2.0) — Diagnostic\n";
        auto metrics_loose = evaluateCEC2015Metrics(candidates, fp, env_raw, 0.05, 2.0);
        std::cout << "    SR       : " << metrics_loose.SR << "%\n"
            << "    ANOF     : " << metrics_loose.ANOF << "  ("
            << static_cast<int>(std::round(metrics_loose.ANOF * n_opt))
            << "/" << n_opt << " found)\n"
            << "    MPR      : " << metrics_loose.MPR << "\n"
            << "    SP       : " << metrics_loose.SP << "\n";

        auto best_sol = solver.getBestConsensusSolution();
        if (!best_sol.empty()) {
            std::cout << ">>> Best consensus solution: [";
            for (size_t i = 0; i < best_sol.size(); i++) {
                std::cout << std::fixed << std::setprecision(4) << best_sol[i];
                if (i + 1 < best_sol.size()) std::cout << ", ";
            }
            std::cout << "]  f = "
                << solver.getBestConsensusFitness() << "\n";
        }

        outputLandscapeSlices(env_raw, vis_base + "/landscape_after");

        std::cout << "\n>>> [MPM-CoEA v2.0] Complete.\n";
    }

} // namespace ofec
#endif // OFEC_MPMCOEA_SOLVER_HPP