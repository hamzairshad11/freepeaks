#ifndef OFEC_MPMCOEA_SOLVER_HPP
#define OFEC_MPMCOEA_SOLVER_HPP

/*
 * BUG FIXES:
 *
 *   F1 – MapXPartyBias clamping (catastrophic)
 *        The original clampToUnitInterval() received function-domain
 *        RESIDUALS in [-200, 200] and forced them to [0, 1].  Every
 *        negative residual (point to the "left" of a peak center)
 *        became 0, giving artificial distance = 0 → fitness ≈ peak_height
 *        everywhere in that half-domain.  This creates a massive flat
 *        plateau that makes DE stagnate immediately.
 *        FIX: Removed all clamping; bias magnitude re-scaled to function
 *        domain (5–20 units instead of 0.03–0.04).
 *        See map_x_partybias.h / map_x_partybias.cpp (separate files).
 *
 *   F2 – Multi-peak-per-subspace → only first peak evaluated
 *        FreePeaks::evaluate_ iterates j ∈ [0, numberObjectives()-1] = [0,0],
 *        so only the first OnePeak in a subspace is ever used.  The original
 *        design placed 7 peaks per subspace, producing 14 "optima" in the
 *        metadata while the landscape had only 2 real peaks (1 per subspace).
 *        This explains ANOF = 1/14 = 0.0714 from the start.
 *        FIX: 8 subspaces, exactly 1 peak each → 8 real optima in landscape.
 *
 *   F3 – Algorithm stagnation (never improves after gen 0)
 *        NBC produced sub-populations of size 5 (at MIN_POP_SIZE).  With
 *        all 5 members in a tight cluster, DE trials are indistinguishable
 *        from the current best.  Context vectors froze after 1 generation.
 *        FIX: larger MIN_POP_SIZE (20), larger pool (8 000), adaptive restart
 *        when no improvement for STAGNATION_LIMIT generations, context
 *        updated from ALL sub-pop peaks (not just global best).
 *
 *   F4 – Solution checkpoints only at every 50th generation
 *        FIX: checkpoints every 5 generations.
 *
 */


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

namespace ofec {

    //  NEAREST-BETTER CLUSTERING  (NBC)

    std::vector<std::vector<int>> nearestBetterClustering(
        const std::vector<std::vector<double>>& population,
        const std::vector<double>& fitness,
        double                                  phi = 2.0)
    {
        const int n = static_cast<int>(population.size());
        if (n == 0) return {};

        struct Edge { int from, to; double weight; };
        std::vector<Edge> edges;
        edges.reserve(static_cast<size_t>(n));

        for (int i = 0; i < n; i++) {
            double min_dist = std::numeric_limits<double>::max();
            int    nb_idx = -1;
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
                std::queue<int>  q;
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


    //  MPMMO_Benchmark
    //  Wraps FreePeaks and manages the variable decomposition across parties.
    class MPMMO_Benchmark {
    private:
        FreePeaks* m_problem;
        Environment* m_env;
        int          m_num_dms;
        std::vector<std::vector<int>> m_party_indices;  // [dm] = global var indices

    public:
        MPMMO_Benchmark(FreePeaks* prob, int num_dms, Environment* env)
            : m_problem(prob), m_env(env), m_num_dms(num_dms)
        {
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

        double evaluateJoint(const std::vector<double>& full_sol,
            Environment* env) const
        {
            std::unique_ptr<SolutionBase> sol(m_problem->createSolution());
            auto& sv = dynamic_cast<Solution<VariableVector<Real>>&>(*sol)
                .variable().vector();
            for (size_t i = 0; i < sv.size() && i < full_sol.size(); i++)
                sv[i] = static_cast<Real>(full_sol[i]);
            sol->evaluate(env, false);
            return static_cast<double>(sol->objective(0));
        }

        // Fitness for party dm given its sub-vector and the other parties' contexts.
        double evaluateDM(int                                      dm_id,
            const std::vector<double>& my_vars,
            const std::vector<std::vector<double>>& context,
            Environment* env) const
        {
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

        int          getNumDMs()        const { return m_num_dms; }
        int          getNumVariables()  const { return static_cast<int>(m_problem->numberVariables()); }
        int          getPartyDim(int p) const { return static_cast<int>(m_party_indices[static_cast<size_t>(p)].size()); }
        const std::vector<int>& getPartyIndices(int p) const { return m_party_indices[static_cast<size_t>(p)]; }
        FreePeaks* getProblem() const { return m_problem; }
        Environment* getEnv()     const { return m_env; }
    };



    //  MPM_CoEA – Multiparty Multimodal Co-Evolutionary Algorithm
    class MPM_CoEA {
    private:
        std::shared_ptr<MPMMO_Benchmark> m_bench;
        std::shared_ptr<Environment>     m_env;

        // Algorithm parameters
        int m_pop_size;        // target individuals per competing sub-pop
        int m_med_pop_size;    // mediating (joint) population size
        int m_K;               // max NBC sub-populations per party

        // Larger minimum so DE always has enough distinct individuals
        static constexpr int MIN_POP_SIZE = 20;
        static constexpr int STAGNATION_LIMIT = 60;   // gens before forced restart
        static constexpr double DE_F = 0.65;
        static constexpr double DE_CR = 0.85;

        // Competing populations  [dm][subpop][ind] = DM-specific sub-vector
        std::vector<std::vector<std::vector<std::vector<double>>>> m_competing_pops;
        std::vector<std::vector<std::vector<double>>>              m_competing_fitness;

        // Mediating (joint) population  [ind] = full joint solution
        std::vector<std::vector<double>> m_mediating_pop;
        std::vector<double>              m_mediating_fitness;

        // MULTI-CONTEXT ARCHIVE
        // Each party maintains up to K diverse context vectors.
        // Sub-population sp of party dm evaluates against 
        //   archive[other][sp % archive_size[other]].
        std::vector<std::vector<std::vector<double>>> m_context_archive; // [dm][idx][vars]
        std::vector<std::vector<double>>              m_context_fitness; // [dm][idx]

        // EXPLICIT JOINT OPTIMA ARCHIVE (radius clearing)
        struct FoundOptimum {
            std::vector<double> sol;
            double fitness;
        };
        std::vector<FoundOptimum> m_found_optima;
        static constexpr double CLEARING_R = 0.15; // ~15% of unit hypercube; merges elongated basin

        // Per-party elite archive (kept for injectConsensusCandidates random mixes)
        std::vector<std::vector<std::vector<double>>> m_elite_archive;   // [dm][k][vars]
        std::vector<std::vector<double>>              m_elite_fitness;   // [dm][k]

        std::string m_output_directory;
        int         m_generation = 0;
        int         m_stagnation_counter = 0;
        double      m_best_prev = -std::numeric_limits<double>::max();

        // Convergence log: (generation, best_joint_fitness)
        std::vector<std::pair<int, double>> m_convergence_log;

    public:
        MPM_CoEA(std::shared_ptr<MPMMO_Benchmark> bench,
            std::shared_ptr<Environment>     env,
            int                              pop_size = 40,
            int                              med_pop = 400,
            int                              K = 10,
            const std::string& output_dir = "visualization/solutions")
            : m_bench(bench), m_env(env),
            m_pop_size(pop_size), m_med_pop_size(med_pop), m_K(K),
            m_output_directory(output_dir)
        {
            std::filesystem::create_directories(output_dir);
        }


        //  initialize()
        void initialize() {
            const int total_dim = m_bench->getNumVariables();
            const int num_dms = m_bench->getNumDMs();

            // Pool must be large enough to cover all peaks
            // Rule of thumb: 1000 samples per expected peak × num_dms
            const int POOL_SIZE = std::max(8000, 1000 * total_dim * num_dms);

            std::mt19937 rng(42u);
            std::uniform_real_distribution<double> unif(0.0, 1.0);

            // Uniform random pool in [0,1]^dim
            std::vector<std::vector<double>> pool(
                static_cast<size_t>(POOL_SIZE),
                std::vector<double>(static_cast<size_t>(total_dim)));
            for (auto& sol : pool)
                for (auto& v : sol) v = unif(rng);

            // Joint fitness
            std::vector<double> joint_fitness(static_cast<size_t>(POOL_SIZE));
            for (int i = 0; i < POOL_SIZE; i++)
                joint_fitness[static_cast<size_t>(i)] =
                m_bench->evaluateJoint(pool[static_cast<size_t>(i)], m_env.get());

            {
                double mx = *std::max_element(joint_fitness.begin(), joint_fitness.end());
                double mn = *std::min_element(joint_fitness.begin(), joint_fitness.end());
                std::cout << "[INFO] Pool joint fitness range: ["
                    << std::fixed << std::setprecision(3)
                    << mn << ", " << mx << "]\n";
            }

            // Neutral-context per-party fitness
            // Neutral context: every party at domain midpoint 0.5
            std::vector<std::vector<double>> neutral_ctx(static_cast<size_t>(num_dms));
            for (int dm = 0; dm < num_dms; dm++)
                neutral_ctx[static_cast<size_t>(dm)].assign(
                    static_cast<size_t>(m_bench->getPartyDim(dm)), 0.5);

            // per_fit[dm][i] = party dm's view of pool individual i
            std::vector<std::vector<double>> per_fit(
                static_cast<size_t>(num_dms),
                std::vector<double>(static_cast<size_t>(POOL_SIZE)));

            for (int dm = 0; dm < num_dms; dm++) {
                const int    dm_dim = m_bench->getPartyDim(dm);
                const auto& idx = m_bench->getPartyIndices(dm);
                for (int i = 0; i < POOL_SIZE; i++) {
                    std::vector<double> my_vars(static_cast<size_t>(dm_dim));
                    for (int d = 0; d < dm_dim; d++)
                        my_vars[static_cast<size_t>(d)] =
                        pool[static_cast<size_t>(i)][static_cast<size_t>(idx[d])];
                    per_fit[static_cast<size_t>(dm)][static_cast<size_t>(i)] =
                        m_bench->evaluateDM(dm, my_vars, neutral_ctx, m_env.get());
                }
            }

            // NBC clustering per party
            m_competing_pops.resize(static_cast<size_t>(num_dms));
            m_competing_fitness.resize(static_cast<size_t>(num_dms));
            m_elite_archive.resize(static_cast<size_t>(num_dms));
            m_elite_fitness.resize(static_cast<size_t>(num_dms));

            for (int dm = 0; dm < num_dms; dm++) {
                const int    dm_dim = m_bench->getPartyDim(dm);
                const auto& idx = m_bench->getPartyIndices(dm);

                // Project pool to party dimension for NBC distances
                std::vector<std::vector<double>> dm_pool(
                    static_cast<size_t>(POOL_SIZE),
                    std::vector<double>(static_cast<size_t>(dm_dim)));
                for (int i = 0; i < POOL_SIZE; i++)
                    for (int d = 0; d < dm_dim; d++)
                        dm_pool[static_cast<size_t>(i)][static_cast<size_t>(d)] =
                        pool[static_cast<size_t>(i)][static_cast<size_t>(idx[d])];

                // Use phi=1.5 for slightly tighter clusters (more sub-pops)
                auto clusters = nearestBetterClustering(
                    dm_pool, per_fit[static_cast<size_t>(dm)], 1.5);

                // Sort by best-member fitness (descending)
                std::sort(clusters.begin(), clusters.end(),
                    [&](const std::vector<int>& a, const std::vector<int>& b) {
                        double ba = -1e30, bb = -1e30;
                        for (int ii : a)
                            ba = std::max(ba, per_fit[dm][static_cast<size_t>(ii)]);
                        for (int ii : b)
                            bb = std::max(bb, per_fit[dm][static_cast<size_t>(ii)]);
                        return ba > bb;
                    });

                // Cap at K sub-pops (prevents thousands of singletons)
                if (static_cast<int>(clusters.size()) > m_K)
                    clusters.resize(static_cast<size_t>(m_K));

                m_competing_pops[static_cast<size_t>(dm)].clear();
                m_competing_fitness[static_cast<size_t>(dm)].clear();

                std::uniform_real_distribution<double> noise(-0.02, 0.02);

                for (auto& cluster : clusters) {
                    std::sort(cluster.begin(), cluster.end(),
                        [&](int a, int b) {
                            return per_fit[dm][static_cast<size_t>(a)]
            > per_fit[dm][static_cast<size_t>(b)];
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
                        sub_fit[static_cast<size_t>(k)] =
                            per_fit[dm][static_cast<size_t>(pi)];
                    }

                    // Pad to MIN_POP_SIZE with small perturbations of best member
                    while (static_cast<int>(sub_pop.size()) < MIN_POP_SIZE) {
                        std::vector<double> perturbed = sub_pop[0];
                        for (auto& v : perturbed) {
                            v += noise(rng);
                            if (v < 0.0) v = 0.0;
                            if (v > 1.0) v = 1.0;
                        }
                        sub_pop.push_back(perturbed);
                        sub_fit.push_back(
                            m_bench->evaluateDM(dm, perturbed, neutral_ctx, m_env.get()));
                    }

                    m_competing_pops[static_cast<size_t>(dm)].push_back(std::move(sub_pop));
                    m_competing_fitness[static_cast<size_t>(dm)].push_back(std::move(sub_fit));
                }

                // Build per-party elite archive (top-K best individuals)
                buildEliteArchive(dm, dm_pool, per_fit[dm]);

                std::cout << "[INFO] Party " << dm << ": "
                    << m_competing_pops[dm].size() << " sub-pop(s), sizes:";
                for (const auto& sp : m_competing_pops[dm])
                    std::cout << " " << sp.size();
                std::cout << "\n";
            }

            // Mediating population
            // Select the most diverse high-fitness individuals from the pool
            std::vector<int> joint_order(static_cast<size_t>(POOL_SIZE));
            std::iota(joint_order.begin(), joint_order.end(), 0);
            std::sort(joint_order.begin(), joint_order.end(),
                [&](int a, int b) {
                    return joint_fitness[static_cast<size_t>(a)]
                    > joint_fitness[static_cast<size_t>(b)];
                });

            m_mediating_pop.clear();
            m_mediating_fitness.clear();
            m_mediating_pop.reserve(static_cast<size_t>(m_med_pop_size));
            m_mediating_fitness.reserve(static_cast<size_t>(m_med_pop_size));

            const double diversity_r = 0.02;  // minimum L2 distance in [0,1]^d
            for (int rank_i : joint_order) {
                if (static_cast<int>(m_mediating_pop.size()) >= m_med_pop_size) break;
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
            // Pad with random individuals if needed
            while (static_cast<int>(m_mediating_pop.size()) < m_med_pop_size) {
                std::vector<double> rand_sol(static_cast<size_t>(total_dim));
                for (auto& v : rand_sol) v = unif(rng);
                m_mediating_pop.push_back(rand_sol);
                m_mediating_fitness.push_back(
                    m_bench->evaluateJoint(rand_sol, m_env.get()));
            }

            // Initialize MULTI-CONTEXT archives from sub-population bests
            const int num_dms_init = m_bench->getNumDMs();
            m_context_archive.resize(num_dms_init);
            m_context_fitness.resize(num_dms_init);

            for (int dm = 0; dm < num_dms_init; ++dm) {
                m_context_archive[dm].clear();
                m_context_fitness[dm].clear();
                for (size_t sp = 0; sp < m_competing_pops[dm].size(); ++sp) {
                    if (m_competing_pops[dm][sp].empty()) continue;
                    m_context_archive[dm].push_back(m_competing_pops[dm][sp][0]);
                    m_context_fitness[dm].push_back(m_competing_fitness[dm][sp][0]);
                }
                pruneContextArchiveDiverse(dm, 0.05);
            }

            // Initialize found-optima archive from mediating population
            updateFoundOptimaArchive();

            m_best_prev = getBestConsensusFitness();
            m_stagnation_counter = 0;
            m_convergence_log.clear();
            m_convergence_log.push_back({ 0, m_best_prev });

            std::cout << "[INFO] Init complete | med_pop=" << m_mediating_pop.size()
                << " | best_joint=" << std::fixed << std::setprecision(4)
                << m_best_prev << "\n";
        }


        //  run_one_generation()
        void run_one_generation() {
            ++m_generation;
            const int num_dms = m_bench->getNumDMs();
            const int total_dim = m_bench->getNumVariables();

            // DE on each party's competing sub-populations (MULTI-CONTEXT)
            for (int dm = 0; dm < num_dms; dm++) {
                const int dm_dim = m_bench->getPartyDim(dm);

                for (size_t sp = 0; sp < m_competing_pops[dm].size(); sp++) {
                    auto& pop = m_competing_pops[dm][sp];
                    auto& fit = m_competing_fitness[dm][sp];
                    const int pop_n = static_cast<int>(pop.size());
                    if (pop_n < 4) continue;

                    // Build per-sub-pop context. Each sp sees a different opponent belief
                    std::vector<std::vector<double>> ctx(num_dms);
                    for (int other = 0; other < num_dms; ++other) {
                        if (other == dm) {
                            ctx[other].assign(dm_dim, 0.5); // dummy, ignored by evaluateDM
                        }
                        else {
                            const auto& arch = m_context_archive[other];
                            if (!arch.empty()) {
                                size_t ctx_idx = sp % arch.size();
                                ctx[other] = arch[ctx_idx];
                            }
                            else {
                                ctx[other].assign(m_bench->getPartyDim(other), 0.5);
                            }
                        }
                    }

                    std::mt19937 rng_de(static_cast<unsigned>(
                        m_generation * 100003u + dm * 1009u + static_cast<unsigned>(sp)));
                    std::uniform_real_distribution<double> unif(0.0, 1.0);
                    std::uniform_int_distribution<int>    pick(0, pop_n - 1);

                    for (int i = 0; i < pop_n; i++) {
                        int r1 = i, r2 = i, r3 = i, att = 0;
                        while (r1 == i && ++att < 2000) r1 = pick(rng_de);
                        att = 0; while ((r2 == i || r2 == r1) && ++att < 2000) r2 = pick(rng_de);
                        att = 0; while ((r3 == i || r3 == r1 || r3 == r2) && ++att < 2000) r3 = pick(rng_de);
                        if (r1 == i || r2 == i || r3 == i) continue;

                        std::vector<double> trial(pop[i]);
                        const int jrand = static_cast<int>(unif(rng_de) * dm_dim);
                        for (int d = 0; d < dm_dim; d++) {
                            if (unif(rng_de) < DE_CR || d == jrand) {
                                double v = pop[r1][d] + DE_F * (pop[r2][d] - pop[r3][d]);
                                if (v < 0.0) v = 0.0;
                                if (v > 1.0) v = 1.0;
                                trial[d] = v;
                            }
                        }

                        const double tf = m_bench->evaluateDM(dm, trial, ctx, m_env.get());
                        if (tf >= fit[i]) {
                            pop[i] = std::move(trial);
                            fit[i] = tf;
                        }
                    }

                    // Sort best-first
                    std::vector<size_t> order(pop_n);
                    std::iota(order.begin(), order.end(), 0);
                    std::sort(order.begin(), order.end(),
                        [&](size_t a, size_t b) { return fit[a] > fit[b]; });
                    std::vector<std::vector<double>> s_pop(pop_n);
                    std::vector<double> s_fit(pop_n);
                    for (int k = 0; k < pop_n; k++) {
                        s_pop[k] = pop[order[k]];
                        s_fit[k] = fit[order[k]];
                    }
                    pop = std::move(s_pop);
                    fit = std::move(s_fit);
                }
            }

            // UPDATE CONTEXT ARCHIVES from sub-population bests
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

                // TOP mediating projections (fitness-based)
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

                // RANDOM mediating projections (diversity injection)
                // Ensures some contexts come from random explorers, not just high-fitness consensus individuals.
                std::mt19937 rng_div(static_cast<unsigned>(m_generation * 44444u + dm));
                std::uniform_int_distribution<int> pick_rand(0, static_cast<int>(m_mediating_pop.size()) - 1);
                for (int rand_i = 0; rand_i < 3; ++rand_i) {
                    int ri = pick_rand(rng_div);
                    std::vector<double> ctx(dm_dim);
                    for (int d = 0; d < dm_dim; ++d) ctx[d] = m_mediating_pop[ri][idx[d]];
                    candidates.push_back(std::move(ctx));
                    cand_fit.push_back(m_mediating_fitness[ri]);
                }

                m_context_archive[dm] = std::move(candidates);
                m_context_fitness[dm] = std::move(cand_fit);
                pruneContextArchiveDiverse(dm, 0.05);

                updateEliteArchive(dm);
            }

            // SPECIES-BASED MEDIATING EVOLUTION
            evolveMediatingSpecies();

            // Update explicit optima archive after mediating evolution
            updateFoundOptimaArchive();
            pruneFoundOptimaArchive();          // Collapse 1080 → ~1-8 entries
            
            if (m_generation % 20 == 0) {
                std::cout << "  [archive fitnesses:";
                for (const auto& opt : m_found_optima) {
                    std::cout << " " << std::fixed << std::setprecision(1) << opt.fitness;
                }
                std::cout << "]\n";
            }

            // Inject consensus candidates into mediating population
            injectConsensusCandidates();

            // Enforce per-basin capacity, eject excess to random exploration
            clearMediatingPopulation();

            // Feed diverse mediating projections back into context archives
            mergeContextsFromMediating();

            // Stagnation detection and adaptive restart
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
                    << cur_best << ") — re-seeding competing populations\n";
                restartCompetingPops();
                m_stagnation_counter = 0;
            }

            // Periodic console report
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

        // Polish each archive entry toward exact peak center
        // Runs a mini-DE of 15 individuals for 40 generations around each
        // found optimum.  This bridges the gap between basin-level finding
        // (CLEARING_R=0.18) and CEC peak-center precision (eps_x=0.05).
        void refineFoundOptima() {
            if (m_found_optima.empty()) return;

            std::mt19937 rng(99999u);
            std::uniform_real_distribution<double> noise(-0.03, 0.03);
            const int total_dim = m_bench->getNumVariables();

            for (auto& opt : m_found_optima) {
                // Mini-population around archive entry
                const int mini_n = 15;
                std::vector<std::vector<double>> mini_pop(mini_n, opt.sol);
                std::vector<double> mini_fit(mini_n, opt.fitness);

                // Perturb to create diversity
                for (int i = 1; i < mini_n; ++i) {
                    for (auto& v : mini_pop[i]) {
                        v += noise(rng);
                        if (v < 0.0) v = 0.0;
                        if (v > 1.0) v = 1.0;
                    }
                    mini_fit[i] = m_bench->evaluateJoint(mini_pop[i], m_env.get());
                }

                // Mini DE for 40 gens
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
                        for (int d = 0; d < total_dim; ++d) {
                            if (unif(rng) < DE_CR || d == jrand) {
                                double v = mini_pop[r1][d] + DE_F * (mini_pop[r2][d] - mini_pop[r3][d]);
                                if (v < 0.0) v = 0.0;
                                if (v > 1.0) v = 1.0;
                                trial[d] = v;
                            }
                        }
                        double tf = m_bench->evaluateJoint(trial, m_env.get());
                        if (tf >= mini_fit[i]) {
                            mini_pop[i] = std::move(trial);
                            mini_fit[i] = tf;
                        }
                    }
                }

                // Replace archive with best refined individual
                auto best_it = std::max_element(mini_fit.begin(), mini_fit.end());
                size_t bi = static_cast<size_t>(std::distance(mini_fit.begin(), best_it));
                opt.sol = mini_pop[bi];
                opt.fitness = mini_fit[bi];
            }
        }


        double getBestConsensusFitness() const {
            if (m_mediating_fitness.empty())
                return -std::numeric_limits<double>::max();
            return *std::max_element(m_mediating_fitness.begin(), m_mediating_fitness.end());
        }

        std::vector<double> getBestConsensusSolution() const {
            if (m_mediating_pop.empty()) return {};
            const auto it = std::max_element(m_mediating_fitness.begin(), m_mediating_fitness.end());
            return m_mediating_pop[static_cast<size_t>(
                std::distance(m_mediating_fitness.begin(), it))];
        }


        // Getter for refined archive
        std::vector<std::vector<double>> getFoundOptimaSolutions() const {
            std::vector<std::vector<double>> result;
            for (const auto& opt : m_found_optima) result.push_back(opt.sol);
            return result;
        }


        const std::vector<std::vector<double>>& getMediatingPopulation() const {
            return m_mediating_pop;
        }


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

        //  OUTPUT PARTY-SPECIFIC LANDSCAPE SLICES
        void outputPartyLandscapes(const std::string& output_dir) const {
            std::filesystem::create_directories(output_dir);
            const int res = 200;
            const int num_dms = m_bench->getNumDMs();

            // Contexts to test: neutral, global-best projection, random diverse
            std::vector<std::vector<std::vector<double>>> contexts;
            std::vector<std::string> ctx_names;

            // Neutral context
            std::vector<std::vector<double>> neutral(num_dms);
            for (int dm = 0; dm < num_dms; ++dm)
                neutral[dm].assign(m_bench->getPartyDim(dm), 0.5);
            contexts.push_back(neutral);
            ctx_names.push_back("neutral");

            // Best consensus context (from mediating pop)
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

            // A random context from each party's archive (diversity)
            if (!m_context_archive.empty()) {
                std::vector<std::vector<double>> rand_ctx(num_dms);
                for (int dm = 0; dm < num_dms; ++dm) {
                    if (!m_context_archive[dm].empty())
                        rand_ctx[dm] = m_context_archive[dm].back(); // arbitrary diverse entry
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

                            // Inject the opponent context
                            for (int other = 0; other < num_dms; ++other) {
                                if (other == dm) continue;
                                const auto& oth_idx = m_bench->getPartyIndices(other);
                                const auto& ctx = contexts[ctx_id][other];
                                for (size_t k = 0; k < oth_idx.size(); ++k)
                                    full_sol[oth_idx[k]] = ctx[k];
                            }

                            double fit = m_bench->evaluateJoint(full_sol, m_env.get());
                            f << full_sol[g1] << " " << full_sol[g2] << " " << fit << "\n";
                        }
                    }
                }
            }
        }

    private:

        // Build elite archive: top-K distinct best individuals for party dm
        void buildEliteArchive(int dm,
            const std::vector<std::vector<double>>& pool,
            const std::vector<double>& pool_fit)
        {
            const int dm_dim = m_bench->getPartyDim(dm);
            int POOL_SIZE = static_cast<int>(pool.size());

            std::vector<int> order(static_cast<size_t>(POOL_SIZE));
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(),
                [&](int a, int b) { return pool_fit[a] > pool_fit[b]; });

            m_elite_archive[dm].clear();
            m_elite_fitness[dm].clear();
            const double min_sep = 0.05;  // minimum distance in [0,1]^dm_dim

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

        // Elite archive
        void updateEliteArchive(int dm) {
            for (size_t sp = 0; sp < m_competing_pops[dm].size(); sp++) {
                if (m_competing_pops[dm][sp].empty()) continue;
                const auto& best_ind = m_competing_pops[dm][sp][0];
                const double best_f = m_competing_fitness[dm][sp][0];

                // Add to archive if better than worst and diverse enough
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

        // Inject all cross-party combination candidates into mediating population
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
            // For 2 parties this yields 3x3 = 9 high-quality joint candidates.
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


            // Random elite mixes for stochastic exploration
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


            // Inject random subspace combinations
            // For each party, take a random sub-pop best and combine with
            // random context from other party. This creates "wild" combinations
            // that may hit unexplored joint optima.
            std::mt19937 rng_wild(static_cast<unsigned>(m_generation * 88888u));
            for (int wild = 0; wild < 20; ++wild) {
                std::vector<double> cand(total_dim, 0.5);
                for (int dm = 0; dm < num_dms; ++dm) {
                    const auto& idx = m_bench->getPartyIndices(dm);
                    if (m_competing_pops[dm].empty() || m_context_archive[dm].empty()) continue;

                    // 50% chance: use sub-pop best, 50%: use random context
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

        void insertIfBetter(const std::vector<double>& cand) {
            const double raw_cf = m_bench->evaluateJoint(cand, m_env.get());
            const double cf = penalizedFitness(cand, raw_cf);
            auto worst_it = std::min_element(m_mediating_fitness.begin(), m_mediating_fitness.end());
            if (cf > *worst_it) {
                const size_t wi = static_cast<size_t>(
                    std::distance(m_mediating_fitness.begin(), worst_it));
                m_mediating_pop[wi] = cand;
                m_mediating_fitness[wi] = raw_cf; // store raw fitness
            }
        }

        // Re-seed competing populations from mediating population
        void restartCompetingPops() {
            const int num_dms = m_bench->getNumDMs();
            std::mt19937 rng_rs(static_cast<unsigned>(m_generation * 55555u));
            std::uniform_real_distribution<double> unif(0.0, 1.0);

            for (int dm = 0; dm < num_dms; dm++) {
                const int dm_dim = m_bench->getPartyDim(dm);
                const auto& idx = m_bench->getPartyIndices(dm);

                // Keep only the best sub-population, discard the rest
                if (!m_competing_pops[dm].empty()) {
                    auto best_sp = std::move(m_competing_pops[dm][0]);
                    auto best_ft = std::move(m_competing_fitness[dm][0]);
                    m_competing_pops[dm].clear();
                    m_competing_fitness[dm].clear();
                    m_competing_pops[dm].push_back(std::move(best_sp));
                    m_competing_fitness[dm].push_back(std::move(best_ft));
                }

                // Re-seed remaining sub-pops with repulsion from found optima
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
                            for (const auto& opt : m_found_optima) {
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

                        // Evaluate against NEUTRAL context (0.5) to give new sub-pops a fair
                        // chance to establish in their own basins before facing the opponent's hardened consensus.
                        std::vector<std::vector<double>> neutral_ctx(num_dms);
                        for (int d = 0; d < num_dms; ++d)
                            neutral_ctx[d].assign(m_bench->getPartyDim(d), 0.5);
                        sub_fit[k] = m_bench->evaluateDM(dm, sub_pop[k], neutral_ctx, m_env.get());
                    }
                    m_competing_pops[dm].push_back(std::move(sub_pop));
                    m_competing_fitness[dm].push_back(std::move(sub_fit));
                }
            }
        }

        void updateContextFromMediating() {
            const int num_dms = m_bench->getNumDMs();
            m_context_archive.resize(num_dms);
            m_context_fitness.resize(num_dms);

            if (m_mediating_pop.empty()) {
                for (int dm = 0; dm < num_dms; ++dm) {
                    m_context_archive[dm].clear();
                    m_context_fitness[dm].clear();
                    m_context_archive[dm].push_back(std::vector<double>(m_bench->getPartyDim(dm), 0.5));
                    m_context_fitness[dm].push_back(-1e30);
                }
                return;
            }

            std::vector<size_t> order(m_mediating_pop.size());
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(),
                [&](size_t a, size_t b) { return m_mediating_fitness[a] > m_mediating_fitness[b]; });

            for (int dm = 0; dm < num_dms; ++dm) {
                const auto& idx = m_bench->getPartyIndices(dm);
                const int dm_dim = m_bench->getPartyDim(dm);

                m_context_archive[dm].clear();
                m_context_fitness[dm].clear();

                for (size_t rank = 0; rank < order.size() && m_context_archive[dm].size() < static_cast<size_t>(m_K); ++rank) {
                    std::vector<double> ctx(dm_dim);
                    for (int d = 0; d < dm_dim; ++d)
                        ctx[d] = m_mediating_pop[order[rank]][idx[d]];
                    m_context_archive[dm].push_back(std::move(ctx));
                    m_context_fitness[dm].push_back(m_mediating_fitness[order[rank]]);
                }
                pruneContextArchiveDiverse(dm, 0.05);
            }
        }


        /*// Context archive maintenance
        void pruneContextArchive(int dm, double min_sep) {
            auto& arch = m_context_archive[dm];
            auto& fit = m_context_fitness[dm];
            if (arch.empty()) return;

            std::vector<size_t> order(arch.size());
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(),
                [&](size_t a, size_t b) { return fit[a] > fit[b]; });

            std::vector<std::vector<double>> new_arch;
            std::vector<double> new_fit;
            new_arch.reserve(arch.size());
            new_fit.reserve(arch.size());

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
        }*/

        void pruneContextArchiveDiverse(int dm, double min_sep) {
            auto& arch = m_context_archive[dm];
            auto& fit = m_context_fitness[dm];
            if (arch.empty()) return;

            // Sort by fitness descending
            std::vector<size_t> order(arch.size());
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(),
                [&](size_t a, size_t b) { return fit[a] > fit[b]; });

            std::vector<std::vector<double>> new_arch;
            std::vector<double> new_fit;

            // Greedy diversity-aware subset selection
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

        /*void updateFoundOptimaArchive() {
            // Greedily add the best individual from each undiscovered basin
            while (m_found_optima.size() < 8) {
                double best_f = -1e30;
                size_t best_i = 0;
                bool found_new = false;

                for (size_t i = 0; i < m_mediating_pop.size(); ++i) {
                    double min_dist = 1e30;
                    for (const auto& opt : m_found_optima) {
                        double d2 = 0.0;
                        for (size_t d = 0; d < m_mediating_pop[i].size(); ++d) {
                            double diff = m_mediating_pop[i][d] - opt.sol[d];
                            d2 += diff * diff;
                        }
                        min_dist = std::min(min_dist, std::sqrt(d2));
                    }
                    // Must be at least CLEARING_R from ALL existing entries
                    if (min_dist >= CLEARING_R && m_mediating_fitness[i] > best_f) {
                        best_f = m_mediating_fitness[i];
                        best_i = i;
                        found_new = true;
                    }
                }
                if (!found_new) break; // No more distinct basins in population
                m_found_optima.push_back({ m_mediating_pop[best_i], best_f });
            }

            // Upgrade existing entries if a better point is found in the same basin
            for (auto& opt : m_found_optima) {
                for (size_t i = 0; i < m_mediating_pop.size(); ++i) {
                    double d2 = 0.0;
                    for (size_t d = 0; d < m_mediating_pop[i].size(); ++d) {
                        double diff = m_mediating_pop[i][d] - opt.sol[d];
                        d2 += diff * diff;
                    }
                    if (std::sqrt(d2) < CLEARING_R && m_mediating_fitness[i] > opt.fitness) {
                        opt.sol = m_mediating_pop[i];
                        opt.fitness = m_mediating_fitness[i];
                    }
                }
            }
        }*/
        void updateFoundOptimaArchive() {
            // Only add if the candidate is a genuine improvement in an empty basin
            for (size_t i = 0; i < m_mediating_pop.size(); ++i) {
                double min_dist = 1e30;
                for (const auto& opt : m_found_optima) {
                    double d2 = 0.0;
                    for (size_t d = 0; d < m_mediating_pop[i].size(); ++d) {
                        double diff = m_mediating_pop[i][d] - opt.sol[d];
                        d2 += diff * diff;
                    }
                    min_dist = std::min(min_dist, std::sqrt(d2));
                }

                // Basin is empty and fitness is "peak-worthy"
                if (min_dist >= CLEARING_R && m_mediating_fitness[i] > 30.0) {
                    // 30.0 is a floor; adjust to your problem's noise level
                    m_found_optima.push_back({ m_mediating_pop[i], m_mediating_fitness[i] });
                }
            }
            pruneFoundOptimaArchive();
        }


        double penalizedFitness(const std::vector<double>& sol, double raw_fitness) const {
            for (const auto& opt : m_found_optima) {
                double d2 = 0.0;
                for (size_t d = 0; d < sol.size(); ++d) {
                    double diff = sol[d] - opt.sol[d];
                    d2 += diff * diff;
                }
                double dist = std::sqrt(d2);
                if (dist < CLEARING_R) {
                    if (raw_fitness > opt.fitness + 1e-3) {
                        return raw_fitness; // strict improvement allowed
                    }
                    else {
                        // Strong dynamic penalty: proportional to fitness range
                        double penalty = 20.0 * (1.0 - dist / CLEARING_R);
                        return raw_fitness - penalty;
                    }
                }
            }
            return raw_fitness;
        }


        // PRUNE found-optima archive: merge close entries, keep best, cap size
        void pruneFoundOptimaArchive() {
            if (m_found_optima.size() <= 1) return;

            std::vector<size_t> order(m_found_optima.size());
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(),
                [&](size_t a, size_t b) { return m_found_optima[a].fitness > m_found_optima[b].fitness; });

            std::vector<FoundOptimum> pruned;
            for (size_t idx : order) {
                bool too_close = false;
                for (const auto& ex : pruned) {
                    double d2 = 0.0;
                    for (size_t d = 0; d < m_found_optima[idx].sol.size(); ++d) {
                        double diff = m_found_optima[idx].sol[d] - ex.sol[d];
                        d2 += diff * diff;
                    }
                    if (std::sqrt(d2) < CLEARING_R) { too_close = true; break; }
                }
                if (!too_close) {
                    pruned.push_back(m_found_optima[idx]);
                    if (pruned.size() >= 8) break;  // Hard cap at expected number of real peaks
                }
            }
            m_found_optima = std::move(pruned);
        }

        //  SPECIES-BASED MEDIATING EVOLUTION
        //  NBC clusters the 500 joint individuals into species (niches), 
        //  then DE runs independently inside each species.  This prevents global-peak collapse.
        void evolveMediatingSpecies() {
            const int total_dim = m_bench->getNumVariables();
            const int med_n = static_cast<int>(m_mediating_pop.size());
            if (med_n < 8) return;  // safety guard

            // Cluster mediating population by NBC
            auto clusters = nearestBetterClustering(m_mediating_pop, m_mediating_fitness, 2.0);

            // Keep only clusters large enough for DE (need 4 individuals)
            std::vector<std::vector<int>> species;
            for (auto& c : clusters) {
                if (static_cast<int>(c.size()) >= 4)
                    species.push_back(c);
            }
            // If everything collapsed to one cluster, split it
            while (species.size() < 2 && !species.empty()) {
                auto& largest = species[0];
                size_t half = largest.size() / 2;
                std::vector<int> second_half(largest.begin() + static_cast<int>(half), largest.end());
                largest.erase(largest.begin() + static_cast<int>(half), largest.end());
                species.push_back(second_half);
            }
            if (species.empty()) return;

            // Evolve each species independently (island model)
            std::mt19937 rng_spec(static_cast<unsigned>(m_generation * 777777u));
            std::uniform_real_distribution<double> unif(0.0, 1.0);

            for (const auto& sp : species) {
                const int sp_size = static_cast<int>(sp.size());
                if (sp_size < 4) continue;

                std::uniform_int_distribution<int> pick(0, sp_size - 1);

                for (int idx = 0; idx < sp_size; ++idx) {
                    int global_i = sp[idx];
                    int r1 = idx, r2 = idx, r3 = idx, att = 0;
                    while (r1 == idx && ++att < 2000) r1 = pick(rng_spec);
                    att = 0; while ((r2 == idx || r2 == r1) && ++att < 2000) r2 = pick(rng_spec);
                    att = 0; while ((r3 == idx || r3 == r1 || r3 == r2) && ++att < 2000) r3 = pick(rng_spec);
                    if (r1 == idx || r2 == idx || r3 == idx) continue;

                    int g1 = sp[r1], g2 = sp[r2], g3 = sp[r3];

                    std::vector<double> trial(m_mediating_pop[global_i]);
                    const int jrand = static_cast<int>(unif(rng_spec) * total_dim);
                    for (int d = 0; d < total_dim; d++) {
                        if (unif(rng_spec) < DE_CR || d == jrand) {
                            double v = m_mediating_pop[g1][d] + DE_F * (m_mediating_pop[g2][d] - m_mediating_pop[g3][d]);
                            v = std::max(0.0, std::min(1.0, v));
                            trial[d] = v;
                        }
                    }

                    const double raw_tf = m_bench->evaluateJoint(trial, m_env.get());
                    const double pen_tf = penalizedFitness(trial, raw_tf);
                    const double pen_cur = penalizedFitness(m_mediating_pop[global_i], m_mediating_fitness[global_i]);

                    if (pen_tf >= pen_cur) {
                        m_mediating_pop[global_i] = std::move(trial);
                        m_mediating_fitness[global_i] = raw_tf;
                    }
                }
            }
        }

        void clearMediatingPopulation() {
            // Only replace the single worst individual per generation
            // This maintains exploration without destroying weak niches
            auto worst_it = std::min_element(m_mediating_fitness.begin(), m_mediating_fitness.end());
            size_t wi = static_cast<size_t>(std::distance(m_mediating_fitness.begin(), worst_it));

            std::mt19937 rng(static_cast<unsigned>(m_generation * 1234567u));
            std::uniform_real_distribution<double> unif(0.0, 1.0);
            for (auto& v : m_mediating_pop[wi]) v = unif(rng);
            m_mediating_fitness[wi] = m_bench->evaluateJoint(m_mediating_pop[wi], m_env.get());
        }


        void mergeContextsFromMediating() {
            const int num_dms = m_bench->getNumDMs();
            for (int dm = 0; dm < num_dms; ++dm) {
                const auto& idx = m_bench->getPartyIndices(dm);
                const int dm_dim = m_bench->getPartyDim(dm);

                // Take top 5 fitness + 5 random + 5 from NON-GLOBAL basins
                std::vector<size_t> order(m_mediating_pop.size());
                std::iota(order.begin(), order.end(), 0);
                std::sort(order.begin(), order.end(),
                    [&](size_t a, size_t b) { return m_mediating_fitness[a] > m_mediating_fitness[b]; });

                // Top 5 by fitness
                for (int r = 0; r < 5 && r < static_cast<int>(order.size()); ++r) {
                    std::vector<double> ctx(dm_dim);
                    for (int d = 0; d < dm_dim; ++d) ctx[d] = m_mediating_pop[order[r]][idx[d]];
                    m_context_archive[dm].push_back(std::move(ctx));
                    m_context_fitness[dm].push_back(m_mediating_fitness[order[r]]);
                }

                // 5 random (exploration)
                std::mt19937 rng_div(static_cast<unsigned>(m_generation * 44444u + dm));
                std::uniform_int_distribution<int> pick_rand(0, static_cast<int>(m_mediating_pop.size()) - 1);
                for (int rand_i = 0; rand_i < 5; ++rand_i) {
                    int ri = pick_rand(rng_div);
                    std::vector<double> ctx(dm_dim);
                    for (int d = 0; d < dm_dim; ++d) ctx[d] = m_mediating_pop[ri][idx[d]];
                    m_context_archive[dm].push_back(std::move(ctx));
                    m_context_fitness[dm].push_back(m_mediating_fitness[ri]);
                }

                // 5 from non-global basins (if any exist)
                if (!m_found_optima.empty()) {
                    auto global_it = std::max_element(m_found_optima.begin(), m_found_optima.end(),
                        [](const FoundOptimum& a, const FoundOptimum& b) { return a.fitness < b.fitness; });
                    for (const auto& opt : m_found_optima) {
                        if (opt.fitness >= global_it->fitness - 1.0) continue; // skip global
                        // Find mediating individual closest to this non-global optimum
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

                pruneContextArchiveDiverse(dm, 0.05);
            }
        }
    };


    //  CEC-2015 Metrics
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
        double       eps_x = 0.05,   // tightened: must be within 5% of domain
        double       eps_f = 2.0)    // within 2 fitness units of true optimum
    {
        CEC2015Metrics m;
        const auto* opt = problem->optima();
        if (!opt || opt->numberSolutions() == 0) return m;

        const int n_opt = static_cast<int>(opt->numberSolutions());
        std::vector<bool> found(static_cast<size_t>(n_opt), false);
        int total_found = 0;

        for (const auto& cand : candidates) {
            std::unique_ptr<SolutionBase> sol(problem->createSolution());
            auto& sv = dynamic_cast<Solution<VariableVector<Real>>&>(*sol)
                .variable().vector();
            for (size_t d = 0; d < cand.size() && d < sv.size(); d++)
                sv[d] = static_cast<Real>(cand[d]);
            sol->evaluate(env, false);
            const double cf = sol->objective(0);

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


    //  createComplex4D2PartyProblem
    //
    //  DESIGN:
    // 
    //  1. dim = 4, num_parties = 2, party 0 → {x0,x2}, party 1 → {x1,x3}.
    //     Each party controls 2 dimensions, requiring genuine cooperation.
    //
    //  2. 8 subspaces (one per peak) → 8 real optima.
    //     With FIX F2, each subspace contains exactly 1 OnePeak.
    //
    //  3. Peak heights and shapes vary to create basin-width diversity:
    //       s1 (linear),  s2 (exp),  s3 (sqrt),  s4 (1/d),
    //       s5 (quadratic), s6 (double-exp), s7 (cosine hat)
    //
    //  4. Per-subspace transforms:
    //       – MapXIllConditioning  : axis stretching (varies per subspace)
    //       – MapXIrregularity     : BBOB-style non-linearity
    //       – MapXAssymetrix       : asymmetric valleys (varies per subspace)
    //
    //  5. Subspace proportions (weights) match the landscape topology:
    //     global optimum gets the largest region, deceptive peaks smaller ones.

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

        // Problem dimensions
        const int numDim = 4;   // 4 variables: 2 per party
        const int numObj = 1;   // single-objective multimodal
        const int numCon = 0;

        std::shared_ptr<ofec::Random> rnd(new Random(0.123));

        FreePeaks::registerFP();
        const std::string fp_name = "free_peaks";
        std::shared_ptr<Environment> env(generateEnvironmentByFactory(fp_name));
        env->recordInputParameters();
        env->initialize();
        env->setProblem(generateProblemByFactory(fp_name));
        auto* freepeak = dynamic_cast<FreePeaks*>(env->problem());

        ParameterMap fp_pm;
        fp_pm["generation_type"] = std::string("assigned");
        fp_pm["dataFile1"] = problem_dir + problem_name + ".txt";
        freepeak->inputParameters().input(fp_pm);
        freepeak->initialize(env.get());
        freepeak->setRandom(rnd);
        freepeak->setSizes(numDim, numObj, numCon);

        // 8-subspace KD-tree
        // Weights proportional to landscape prominence.
        // The global optimum (peak1) gets the largest region.
        freepeak->setKDtree({
            { "root", {
                { "peak1", 0.20 },   // global optimum – large basin
                { "peak2", 0.16 },   // strong local
                { "peak3", 0.13 },   // prominent local
                { "peak4", 0.12 },   // moderate local
                { "peak5", 0.11 },   // moderate local
                { "peak6", 0.10 },   // shallow local
                { "peak7", 0.10 },   // shallow local
                { "peak8", 0.08 }    // deceptive – narrow but tall
            }}
            });

        // Peak specifications
        // center[] : in function-domain [-100, 100]^4
        //            FunctionBase maps subspace [lo,hi] → [-100,100] per dim,
        //            so (0,0,0,0) is always reachable in any subspace.
        // height   : peak amplitude
        // shape    : OnePeak shape class (s1..s7)
        // condition: ill-conditioning ratio for this subspace (axis stretch)
        // use_asym : toggle MapXAssymetrix (adds asymmetric valleys)
        struct PeakSpec {
            std::vector<double> center;
            double   height;
            std::string shape;
            double   condition;
            bool     use_asym;
            int         bias_party;  // 0 or 1 — which party's bias distorts this subspace
            double      magnitude;   // bias shift amplitude in function-domain units
            double      rotation;    // Givens rotation angle (radians)
            double      pb_cond;     // condition inside MapXPartyBias
        };

        
        const std::vector<PeakSpec> peaks = {
            // Global optimum: symmetric, no party bias
            { {  0,   0,   0,   0  }, 90.0, "s1",   1.0,  false, 0,  0.0, 0.0,   1.0 },
            // Strong local: Party 0 bias, mild shift+rotation
            { { 40,  -30,  35, -25  }, 75.0, "s2",  10.0,  true,  0, 15.0, 0.4,  10.0 },
            // Prominent local: Party 1 bias
            { {-50,   40, -45,  35  }, 65.0, "s7",  10.0,  false, 1, 12.0, 0.6,  10.0 },
            // Cooperative peak: moderate stretch, Party 0 bias
            { { 25,   25,  25,  25  }, 55.0, "s3",  50.0,  true,  0, 20.0, 0.3,  50.0 },
            // Quadratic basin: wide, Party 1 bias (mild)
            { {-30,  -30, -30, -30  }, 50.0, "s5",   1.0,  false, 1,  8.0, 0.0,   1.0 },
            // Exponential decay: sharp, Party 1 bias
            { { 60,  -50,  55, -45  }, 42.0, "s2", 100.0,  true,  1, 18.0, 0.5, 100.0 },
            // Flat-plateau: broad, Party 0 bias
            { {-20,   60, -15,  55  }, 35.0, "s4",   1.0,  false, 0,  5.0, 0.0,   1.0 },
            // Deceptive: narrow, Party 1 bias, strong distortion
            { { 80,   80,  75,  80  }, 80.0, "s7", 500.0,  true,  1, 25.0, 0.8, 500.0 }
        };

        // Per-subspace ill-conditioning strengths — must match EnumeratedReal set.
        // These are redundant with PeakSpec::condition above but kept for clarity.
        const std::vector<double> conditions = { 1.0, 10.0, 10.0, 50.0, 1.0, 100.0, 1.0, 500.0 };

        const std::vector<std::string> sp_names = {
            "peak1","peak2","peak3","peak4","peak5","peak6","peak7","peak8"
        };

        // Build each subproblem
        for (size_t s = 0; s < sp_names.size(); s++) {
            const std::string& sname = sp_names[s];
            const PeakSpec& pk = peaks[s];

            ParameterMap spm;
            spm["subspace"] = sname;
            spm["generation_type"] = std::string("assigned");
            spm["dataFile1"] = sub_dir + problem_name + "_" + sname + ".txt";

            auto subpro(Subproblem::create());
            subpro->initialize(spm, freepeak);

            // Distance: Euclidean
            {
                auto dis(FactoryFP<DistanceBase>::produce("Euclidean"));
                ParameterMap dp;
                dis->initialize(freepeak, sname, dp);
                subpro->setDistance(dis);
            }

            // Function: one_peak with exactly 1 peak
            {
                ParameterMap fun_p;
                fun_p["generation_type"] = std::string("assigned");
                fun_p["dataFile1"] = sub_dir + problem_name + "_" + sname + "_fn.txt";

                auto func(FactoryFP<FunctionBase>::produce("one_peak"));
                func->initialize(freepeak, sname, fun_p);
                auto* opf = dynamic_cast<OnePeakFunction*>(func);

                // ONE OnePeak per subspace
                auto onepeak(FactoryFP<OnePeakBase>::produce(pk.shape));
                ParameterMap op;
                op["center_type"] = std::string("assigned");
                op["height"] = static_cast<Real>(pk.height);
                std::vector<Real> center_r(pk.center.size());
                for (size_t d = 0; d < pk.center.size(); d++) {
                    // Clamp to valid function domain [-100, 100]
                    center_r[d] = static_cast<Real>(
                        std::max(-100.0, std::min(100.0, pk.center[d])));
                }
                op["center_postion"] = center_r;
                onepeak->initialize(freepeak, sname, op);
                opf->addOnePeaks(onepeak);

                subpro->setFunction(func);
            }

            // X transforms: full landscape differentiation chain
            //
            // TRANSFORM ORDER (applied left-to-right, innermost first):
            //   1. MapXPartyBias      – party-specific coordinate distortion
            //                           THIS IS THE KEY MPMMO TRANSFORM.
            //                           Different parties perceive the same joint
            //                           solution from rotated/shifted perspectives,
            //                           creating genuine cooperation difficulty.
            //   2. MapXIllConditioning – axis stretching (subspace-level)
            //   3. MapXIrregularity   – BBOB-style non-smooth oscillation
            //   4. MapXAssymetrix     – asymmetric valleys (selected subspaces)
            //
            // This 4-layer chain produces landscapes where:
            //   • Parties have asymmetric local optima (via MapXPartyBias)
            //   • Some dimensions are harder to search than others (ill-cond)
            //   • The landscape has irregular, non-smooth structure (irregularity)
            //   • Some peaks have steep sides on one face only (asymmetry)
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
                // Must be an exact member of the EnumeratedReal set:
                // {1, 10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000}
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
                // Slightly different asymmetry strength per subspace
                tp["alpha"] = static_cast<Real>(0.2 + 0.05 * static_cast<double>(s));
                t->initialize(freepeak, sname, tp);
                subpro->addVariableTransform(t);
            }

            // Deceptive oscillation: creates false local attractors
            {
                auto t(FactoryFP<X_TransformBase>::produce("MapXDeceptive"));
                ParameterMap tp;
                tp["alpha"] = static_cast<Real>(0.05 + 0.02 * static_cast<double>(s));
                t->initialize(freepeak, sname, tp);
                subpro->addVariableTransform(t);
            }

            // Cross-dimensional linkage: non-separable coupling between dimensions
            {
                auto t(FactoryFP<X_TransformBase>::produce("MapXLinkage"));
                ParameterMap tp;
                tp["beta"] = static_cast<Real>(0.1);
                t->initialize(freepeak, sname, tp);
                subpro->addVariableTransform(t);
            }

            // Y transform: objective normalisation
            {
                auto ot(FactoryFP<Y_TransformBase>::produce("map_objective"));
                ParameterMap otp;
                ot->initialize(freepeak, sname, otp);
                subpro->addObjectiveTransform(ot);
            }

            freepeak->setSubproblem(sname, subpro);
        }

        // Bind data and persist to files
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
        // Optima count is printed in runMPMCoEAExperiment() after the problem
        // is cleanly reloaded from the just-written files.
    }


    //  outputLandscapeSlices  (2-D projections for visualization)
    void outputLandscapeSlices(Environment* env, const std::string& output_dir) {
        std::filesystem::create_directories(output_dir);
        const int res = 200;

        auto* con_pro = CAST_CONOP(env->problem());
        const int dim = static_cast<int>(con_pro->numberVariables());

        std::cout << ">>> Sampling landscape slices (res=" << res << ")..." << std::flush;

        // For each pair of dimensions, fix all others at 0.5
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
                        auto& v = dynamic_cast<Solution<VariableVector<Real>>&>(*sol)
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

    //  runMPMCoEAExperiment
    void runMPMCoEAExperiment() {
        ofec::registerInstance();

        ofec::g_working_directory =
            R"(E:\HITSZ\Research\Multimodal_Multiparty_Optimization\ThesisProject\Data\ofec_data_new\)";
        for (auto& c : ofec::g_working_directory)
            if (c == '\\') c = '/';

        const std::string vis_base =
            R"(E:\HITSZ\Research\Multimodal_Multiparty_Optimization\ThesisProject\Visualization)";

        // Create benchmark
        const std::string problem_name = "complex_4d_2p";
        std::cout << "\n>>> [MPM-CoEA] Creating 4D 2-party benchmark...\n";
        createComplex4D2PartyProblem(problem_name);

        // Load problem
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

        // Pre-run landscape visualisation
        outputLandscapeSlices(env_raw, vis_base + "/landscape_before");

        // Build benchmark and solver
        auto bench = std::make_shared<MPMMO_Benchmark>(fp, 2, env_raw);
        auto env_sp = std::shared_ptr<Environment>(env_raw, [](Environment*) {});

        // Parameters tuned for 4D, 2 parties, 8 peaks:
        //   pop_size   = 40  (individuals per competing sub-pop)
        //   med_pop    = 500 (joint population – ~60 per peak)
        //   K          = 10  (max sub-populations per party)
        MPM_CoEA solver(bench, env_sp,
            /*pop_size*/ 40,
            /*med_pop*/  500,
            /*K*/        10,
            vis_base + "/solutions");

        solver.initialize();

        // VISUALIZE PARTY-SPECIFIC LANDSCAPES (before evolution)
        solver.outputPartyLandscapes(vis_base + "/party_landscapes_before");

        // Save gen-0 solutions
        solver.saveSolutionsAtGeneration(0);
        std::cout << "Gen   0 | " << std::fixed << std::setprecision(4)
            << solver.getBestConsensusFitness() << "\n";

        // Evolution loop
        const int MAX_GEN = 500;

        for (int gen = 1; gen <= MAX_GEN; gen++) {
            solver.run_one_generation();

            // Console log every 5 generations
            if (gen % 5 == 0 || gen == MAX_GEN) {
                std::cout << "Gen " << std::setw(4) << gen << " | "
                    << std::fixed << std::setprecision(4)
                    << solver.getBestConsensusFitness() << "\n";
            }

            if (gen % 5 == 0 || gen == MAX_GEN)
                solver.saveSolutionsAtGeneration(gen);
        }

        // VISUALIZE PARTY-SPECIFIC LANDSCAPES (after evolution)
        solver.outputPartyLandscapes(vis_base + "/party_landscapes_after");

        // Save convergence curve
        solver.saveConvergenceCurve();

        // CEC-2015 Metrics
        // Refine archive entries to exact peak centers
        solver.refineFoundOptima();

        // Combine final mediating pop + refined archive for evaluation
        auto candidates = solver.getMediatingPopulation();
        auto refined = solver.getFoundOptimaSolutions();
        candidates.insert(candidates.end(), refined.begin(), refined.end());

        auto metrics = evaluateCEC2015Metrics(
            candidates, fp, env_raw,
            /*eps_x*/ 0.05,
            /*eps_f*/ 2.0);

        std::cout << "\n>>> CEC-2015 Metrics\n"
            << "    n_optima : " << n_opt << "\n"
            << "    SR       : " << metrics.SR << "%\n"
            << "    ANOF     : " << metrics.ANOF << "  ("
            << static_cast<int>(std::round(metrics.ANOF * n_opt))
            << "/" << n_opt << " found)\n"
            << "    MPR      : " << metrics.MPR << "\n"
            << "    SP       : " << metrics.SP << "\n";

        // Best solution printout
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

        // Post-run landscape visualisation
        outputLandscapeSlices(env_raw, vis_base + "/landscape_after");

        std::cout << "\n>>> [MPM-CoEA] Complete.\n";
    }

} // namespace ofec
#endif // OFEC_MPMCOEA_SOLVER_HPP