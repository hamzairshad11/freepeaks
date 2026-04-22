#ifndef OFEC_MPMCOEA_SOLVER_HPP
#define OFEC_MPMCOEA_SOLVER_HPP

/*
 * MPM-CoEA: Multiparty Multimodal Co-Evolutionary Algorithm
 *
 * Design principles:
 *   - Each party p maintains up to K competing sub-populations.
 *     Each sub-pop stores ONLY party p's own variable sub-vector.
 *     Fitness = f_p(x_p | context_{-p}).
 *
 *   - A shared mediating population stores FULL joint solutions.
 *     Fitness = TRUE joint objective f([x_0, x_1, ...]).
 *
 *   - Context vectors: one per party, holding best-known values for
 *     that party's variables.  Seeded from the mediating population.
 *
 * Bugs fixed vs. previous version:
 *   B1: Peak centers were in [0,1] space but FreePeaks treats them as
 *       [-100,100] function-domain coordinates.  All peaks appeared at
 *       [-100,-100] causing near-zero fitness everywhere.
 *       Fix: peaks are now specified directly in [-100,100] space.
 *
 *   B2: NBC produced up to 200 singleton sub-pops from a pool of 200.
 *       DE selection loop "while (r1 == i)" with pop_n=1 is infinite.
 *       Fix: NBC clusters are capped at K; each padded to MIN_POP_SIZE.
 *
 *   B3: Mediating fitness reported as 0 at gen-0 because the init
 *       printed stats before the fitness array was fully populated.
 *       Fix: print after the fitness array is filled.
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

namespace ofec {

    //  NEAREST-BETTER CLUSTERING (NBC)
    //  Returns clusters as indices into the population vector.
    std::vector<std::vector<int>> nearestBetterClustering(
        const std::vector<std::vector<double>>& population,
        const std::vector<double>& fitness,
        double phi = 2.0)
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
            // All solutions have equal fitness, each is its own cluster
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
    class MPMMO_Benchmark {
    private:
        FreePeaks* m_problem;
        Environment* m_env;
        int          m_num_dms;
        std::vector<std::vector<int>> m_party_indices; // [dm] = global var indices

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

            std::cout << "[INFO] Variable decomposition:" << std::endl;
            for (int dm = 0; dm < m_num_dms; dm++) {
                std::cout << "  Party " << dm << " -> variables [";
                const auto& idx = m_party_indices[static_cast<size_t>(dm)];
                for (size_t i = 0; i < idx.size(); i++) {
                    std::cout << idx[i];
                    if (i + 1 < idx.size()) std::cout << ", ";
                }
                std::cout << "]" << std::endl;
            }
        }

        //  evaluateJoint: TRUE shared objective at a full joint solution.
        double evaluateJoint(const std::vector<double>& full_sol,
            Environment* env) const
        {
            std::unique_ptr<SolutionBase> sol(m_problem->createSolution());
            auto& sol_vec =
                dynamic_cast<Solution<VariableVector<Real>>&>(*sol).variable().vector();
            for (size_t i = 0; i < sol_vec.size() && i < full_sol.size(); i++)
                sol_vec[i] = static_cast<Real>(full_sol[i]);
            sol->evaluate(env, false);
            return static_cast<double>(sol->objective(0));
        }

        //  evaluateDM: party dm's view of fitness.
        //    my_vars  : length = party_dim[dm]  (DM-specific sub-vector)
        //    context  : context[p] = best-known values for party p's variables
        //               (length = party_dim[p] each)
        double evaluateDM(int                                     dm_id,
            const std::vector<double>& my_vars,
            const std::vector<std::vector<double>>& context,
            Environment* env) const
        {
            const int total_dim = static_cast<int>(m_problem->numberVariables());
            std::vector<double> full_sol(static_cast<size_t>(total_dim), 0.5);

            // This DM's contribution
            const auto& my_idx = m_party_indices[static_cast<size_t>(dm_id)];
            for (size_t i = 0; i < my_idx.size() && i < my_vars.size(); i++)
                full_sol[static_cast<size_t>(my_idx[i])] = my_vars[i];

            // Other parties from context
            for (int other = 0; other < m_num_dms; other++) {
                if (other == dm_id) continue;
                const auto& oth_idx = m_party_indices[static_cast<size_t>(other)];
                const auto& ctx = context[static_cast<size_t>(other)];
                for (size_t i = 0; i < oth_idx.size() && i < ctx.size(); i++)
                    full_sol[static_cast<size_t>(oth_idx[i])] = ctx[i];
            }

            return evaluateJoint(full_sol, env);
        }

        int          getNumDMs()       const { return m_num_dms; }
        int          getNumVariables() const {
            return static_cast<int>(m_problem->numberVariables());
        }
        int          getPartyDim(int p) const {
            return static_cast<int>(m_party_indices[static_cast<size_t>(p)].size());
        }
        const std::vector<int>& getPartyIndices(int p) const {
            return m_party_indices[static_cast<size_t>(p)];
        }
        FreePeaks* getProblem() const { return m_problem; }
        Environment* getEnv()     const { return m_env; }
    };

    //  MPM_CoEA
    class MPM_CoEA {
    private:
        std::shared_ptr<MPMMO_Benchmark> m_bench;
        std::shared_ptr<Environment>     m_env;

        // Algorithm parameters
        int m_pop_size;       // target individuals per competing sub-pop
        int m_med_pop_size;   // mediating population size
        int m_K;              // max NBC sub-populations per party

        static constexpr int MIN_POP_SIZE = 5; // minimum sub-pop size for DE

        // Competing populations [dm][subpop][ind] = DM-specific sub-vector
        std::vector<std::vector<std::vector<std::vector<double>>>> m_competing_pops;
        std::vector<std::vector<std::vector<double>>>              m_competing_fitness;

        // Mediating population [ind] = full joint solution
        std::vector<std::vector<double>> m_mediating_pop;
        std::vector<double>              m_mediating_fitness;

        // Context vectors [dm] = best-known DM-specific sub-vector
        std::vector<std::vector<double>> m_context_vectors;

        std::string m_output_directory;
        int         m_generation = 0;

    public:
        MPM_CoEA(std::shared_ptr<MPMMO_Benchmark> bench,
            std::shared_ptr<Environment>     env,
            int                              pop_size = 50,
            int                              med_pop = 200,
            int                              K = 5,
            const std::string& output_dir = "visualization/solutions")
            : m_bench(bench), m_env(env),
            m_pop_size(pop_size), m_med_pop_size(med_pop), m_K(K),
            m_output_directory(output_dir)
        {
            std::filesystem::create_directories(output_dir);
        }

        //  initialize
        void initialize() {
            const int total_dim = m_bench->getNumVariables();
            const int num_dms = m_bench->getNumDMs();

            // Pool size: large enough for NBC to see the landscape structure
            const int POOL_SIZE = std::max(500, 50 * total_dim * num_dms);

            std::mt19937 rng(42u);
            std::uniform_real_distribution<double> unif(0.0, 1.0);

            // 1. Random pool in [0,1]^dim
            std::vector<std::vector<double>> pool(
                static_cast<size_t>(POOL_SIZE),
                std::vector<double>(static_cast<size_t>(total_dim)));
            for (auto& sol : pool)
                for (auto& v : sol) v = unif(rng);

            // 2. Joint fitness of every pool individual
            std::vector<double> joint_fitness(static_cast<size_t>(POOL_SIZE));
            for (int i = 0; i < POOL_SIZE; i++)
                joint_fitness[static_cast<size_t>(i)] =
                m_bench->evaluateJoint(pool[static_cast<size_t>(i)], m_env.get());

            {
                double max_f = *std::max_element(joint_fitness.begin(), joint_fitness.end());
                double min_f = *std::min_element(joint_fitness.begin(), joint_fitness.end());
                std::cout << "[INFO] Pool joint fitness range: ["
                    << std::fixed << std::setprecision(3)
                    << min_f << ", " << max_f << "]" << std::endl;
            }

            // 3. Per-party fitness (neutral context: all others at 0.5)
            std::vector<std::vector<double>> neutral_ctx(
                static_cast<size_t>(num_dms));
            for (int dm = 0; dm < num_dms; dm++)
                neutral_ctx[static_cast<size_t>(dm)].assign(
                    static_cast<size_t>(m_bench->getPartyDim(dm)), 0.5);

            // per_fit[dm][i] = party-dm's view of pool individual i
            std::vector<std::vector<double>> per_fit(
                static_cast<size_t>(num_dms),
                std::vector<double>(static_cast<size_t>(POOL_SIZE)));

            for (int dm = 0; dm < num_dms; dm++) {
                const int dm_dim = m_bench->getPartyDim(dm);
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

            // 4. NBC clustering per party; cap at K sub-pops each
            m_competing_pops.resize(static_cast<size_t>(num_dms));
            m_competing_fitness.resize(static_cast<size_t>(num_dms));

            for (int dm = 0; dm < num_dms; dm++) {
                const int    dm_dim = m_bench->getPartyDim(dm);
                const auto& idx = m_bench->getPartyIndices(dm);

                // Build party-specific pool projection (for NBC distances)
                std::vector<std::vector<double>> dm_pool(
                    static_cast<size_t>(POOL_SIZE),
                    std::vector<double>(static_cast<size_t>(dm_dim)));
                for (int i = 0; i < POOL_SIZE; i++)
                    for (int d = 0; d < dm_dim; d++)
                        dm_pool[static_cast<size_t>(i)][static_cast<size_t>(d)] =
                        pool[static_cast<size_t>(i)][static_cast<size_t>(idx[d])];

                auto clusters = nearestBetterClustering(
                    dm_pool, per_fit[static_cast<size_t>(dm)], 2.0);

                // Sort clusters by best member fitness (descending)
                std::sort(clusters.begin(), clusters.end(),
                    [&](const std::vector<int>& a, const std::vector<int>& b) {
                        double fa = *std::max_element(a.begin(), a.end(),
                            [&](int x, int y) {
                                return per_fit[static_cast<size_t>(dm)][static_cast<size_t>(x)]
                                    < per_fit[static_cast<size_t>(dm)][static_cast<size_t>(y)];
                            });
                        double fb = *std::max_element(b.begin(), b.end(),
                            [&](int x, int y) {
                                return per_fit[static_cast<size_t>(dm)][static_cast<size_t>(x)]
                                    < per_fit[static_cast<size_t>(dm)][static_cast<size_t>(y)];
                            });
                        // Use per_fit directly for the cluster best
                        double best_a = per_fit[static_cast<size_t>(dm)][static_cast<size_t>(a[0])];
                        for (int ii : a)
                            best_a = std::max(best_a, per_fit[static_cast<size_t>(dm)][static_cast<size_t>(ii)]);
                        double best_b = per_fit[static_cast<size_t>(dm)][static_cast<size_t>(b[0])];
                        for (int ii : b)
                            best_b = std::max(best_b, per_fit[static_cast<size_t>(dm)][static_cast<size_t>(ii)]);
                        return best_a > best_b;
                    });

                // Keep top-K clusters only (prevents thousands of sub-pops)
                if (static_cast<int>(clusters.size()) > m_K)
                    clusters.resize(static_cast<size_t>(m_K));

                m_competing_pops[static_cast<size_t>(dm)].clear();
                m_competing_fitness[static_cast<size_t>(dm)].clear();

                for (auto& cluster : clusters) {
                    // Sort cluster by fitness (descending)
                    std::sort(cluster.begin(), cluster.end(),
                        [&](int a, int b) {
                            return per_fit[static_cast<size_t>(dm)][static_cast<size_t>(a)]
                    > per_fit[static_cast<size_t>(dm)][static_cast<size_t>(b)];
                        });

                    // Take up to m_pop_size but at least MIN_POP_SIZE
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
                            per_fit[static_cast<size_t>(dm)][static_cast<size_t>(pi)];
                    }

                    std::uniform_real_distribution<double> noise(-0.05, 0.05);
                    while (static_cast<int>(sub_pop.size()) < MIN_POP_SIZE) {
                        std::vector<double> perturbed = sub_pop[0];
                        for (auto& v : perturbed) {
                            v += noise(rng);
                            if (v < 0.0) v = 0.0;
                            if (v > 1.0) v = 1.0;
                        }
                        sub_pop.push_back(perturbed);
                        // Re-evaluate the perturbed individual
                        sub_fit.push_back(
                            m_bench->evaluateDM(dm, perturbed, neutral_ctx, m_env.get()));
                    }

                    m_competing_pops[static_cast<size_t>(dm)].push_back(std::move(sub_pop));
                    m_competing_fitness[static_cast<size_t>(dm)].push_back(std::move(sub_fit));
                }

                std::cout << "[INFO] Party " << dm << ": "
                    << m_competing_pops[static_cast<size_t>(dm)].size()
                    << " sub-pop(s), sizes: ";
                for (const auto& sp : m_competing_pops[static_cast<size_t>(dm)])
                    std::cout << sp.size() << " ";
                std::cout << std::endl;
            }

            // 5. Mediating population (full joint, highest diverse fitness)
            // Sort pool by joint fitness descending, filter by diversity
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

            const double diversity_r = 0.03; // minimum L2 distance in [0,1]^d
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
                    m_mediating_fitness.push_back(
                        joint_fitness[static_cast<size_t>(rank_i)]);
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

            // 6. Context vectors from best mediating individual
            updateContextFromMediating();

            std::cout << "[INFO] Init complete | med_pop=" << m_mediating_pop.size()
                << " | best_joint=" << std::fixed << std::setprecision(4)
                << getBestConsensusFitness() << std::endl;
        }

        //  run_one_generation
        void run_one_generation() {
            ++m_generation;
            const int num_dms = m_bench->getNumDMs();
            const int total_dim = m_bench->getNumVariables();

            //  Step 1: DE on each party's competing sub-populations.
            //  Operates ONLY on that party's own sub-vector.
            for (int dm = 0; dm < num_dms; dm++) {
                const int dm_dim = m_bench->getPartyDim(dm);

                for (size_t sp = 0;
                    sp < m_competing_pops[static_cast<size_t>(dm)].size(); sp++)
                {
                    auto& pop = m_competing_pops[static_cast<size_t>(dm)][sp];
                    auto& fit = m_competing_fitness[static_cast<size_t>(dm)][sp];
                    const int pop_n = static_cast<int>(pop.size());

                    // Safety guard: DE requires pop_n >= 4 (FIX B2)
                    if (pop_n < 4) continue;

                    std::mt19937 rng_de(static_cast<unsigned>(
                        m_generation * 100003 + dm * 1009 + static_cast<int>(sp)));
                    std::uniform_real_distribution<double> unif(0.0, 1.0);
                    std::uniform_int_distribution<int>    pick(0, pop_n - 1);

                    const double F = 0.5;
                    const double CR = 0.9;

                    for (int i = 0; i < pop_n; i++) {
                        // Three distinct mutually-different indices != i
                        int r1 = i, r2 = i, r3 = i;
                        int attempts = 0;
                        while (r1 == i && ++attempts < 1000) r1 = pick(rng_de);
                        attempts = 0;
                        while ((r2 == i || r2 == r1) && ++attempts < 1000) r2 = pick(rng_de);
                        attempts = 0;
                        while ((r3 == i || r3 == r1 || r3 == r2) && ++attempts < 1000) r3 = pick(rng_de);
                        if (r1 == i || r2 == i || r3 == i) continue; // can skip if can't find distinct

                        std::vector<double> trial(pop[static_cast<size_t>(i)]);
                        const int jrand = static_cast<int>(unif(rng_de) * dm_dim);

                        for (int d = 0; d < dm_dim; d++) {
                            if (unif(rng_de) < CR || d == jrand) {
                                double v =
                                    pop[static_cast<size_t>(r1)][static_cast<size_t>(d)]
                                    + F * (pop[static_cast<size_t>(r2)][static_cast<size_t>(d)]
                                        - pop[static_cast<size_t>(r3)][static_cast<size_t>(d)]);
                                if (v < 0.0) v = 0.0;
                                if (v > 1.0) v = 1.0;
                                trial[static_cast<size_t>(d)] = v;
                            }
                        }

                        const double tf =
                            m_bench->evaluateDM(dm, trial, m_context_vectors, m_env.get());
                        if (tf >= fit[static_cast<size_t>(i)]) {
                            pop[static_cast<size_t>(i)] = std::move(trial);
                            fit[static_cast<size_t>(i)] = tf;
                        }
                    }

                    // Keep sub-pop sorted: best first
                    std::vector<size_t> order(static_cast<size_t>(pop_n));
                    std::iota(order.begin(), order.end(), 0u);
                    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
                        return fit[a] > fit[b];
                        });
                    std::vector<std::vector<double>> s_pop(static_cast<size_t>(pop_n));
                    std::vector<double>              s_fit(static_cast<size_t>(pop_n));
                    for (int k = 0; k < pop_n; k++) {
                        s_pop[static_cast<size_t>(k)] = pop[order[static_cast<size_t>(k)]];
                        s_fit[static_cast<size_t>(k)] = fit[order[static_cast<size_t>(k)]];
                    }
                    pop = std::move(s_pop);
                    fit = std::move(s_fit);
                }
            }

            //  Step 2: Update context vectors from best competing-pop peaks.
            for (int dm = 0; dm < num_dms; dm++) {
                double              best_fit = -std::numeric_limits<double>::max();
                std::vector<double> best_vars(
                    static_cast<size_t>(m_bench->getPartyDim(dm)), 0.5);

                for (size_t sp = 0;
                    sp < m_competing_pops[static_cast<size_t>(dm)].size(); sp++)
                {
                    if (m_competing_pops[static_cast<size_t>(dm)][sp].empty()) continue;
                    const double f =
                        m_competing_fitness[static_cast<size_t>(dm)][sp][0];
                    if (f > best_fit) {
                        best_fit = f;
                        best_vars = m_competing_pops[static_cast<size_t>(dm)][sp][0];
                    }
                }
                m_context_vectors[static_cast<size_t>(dm)] = best_vars;
            }

            //  Step 3: DE/rand/1/bin on mediating population (TRUE joint f).
            const int med_n = static_cast<int>(m_mediating_pop.size());
            if (med_n < 4) return; // safety guard

            std::mt19937 rng_med(static_cast<unsigned>(m_generation * 999983u));
            std::uniform_real_distribution<double> unif_med(0.0, 1.0);
            std::uniform_int_distribution<int>    pick_med(0, med_n - 1);

            const double F_med = 0.5;
            const double CR_med = 0.9;

            for (int i = 0; i < med_n; i++) {
                int r1 = i, r2 = i, r3 = i;
                int att = 0;
                while (r1 == i && ++att < 1000)            r1 = pick_med(rng_med);
                att = 0;
                while ((r2 == i || r2 == r1) && ++att < 1000)    r2 = pick_med(rng_med);
                att = 0;
                while ((r3 == i || r3 == r1 || r3 == r2) && ++att < 1000) r3 = pick_med(rng_med);
                if (r1 == i || r2 == i || r3 == i) continue;

                std::vector<double> trial(m_mediating_pop[static_cast<size_t>(i)]);
                const int jrand = static_cast<int>(unif_med(rng_med) * total_dim);

                for (int d = 0; d < total_dim; d++) {
                    if (unif_med(rng_med) < CR_med || d == jrand) {
                        double v =
                            m_mediating_pop[static_cast<size_t>(r1)][static_cast<size_t>(d)]
                            + F_med * (m_mediating_pop[static_cast<size_t>(r2)][static_cast<size_t>(d)]
                                - m_mediating_pop[static_cast<size_t>(r3)][static_cast<size_t>(d)]);
                        if (v < 0.0) v = 0.0;
                        if (v > 1.0) v = 1.0;
                        trial[static_cast<size_t>(d)] = v;
                    }
                }

                const double tf = m_bench->evaluateJoint(trial, m_env.get());
                if (tf >= m_mediating_fitness[static_cast<size_t>(i)]) {
                    m_mediating_pop[static_cast<size_t>(i)] = std::move(trial);
                    m_mediating_fitness[static_cast<size_t>(i)] = tf;
                }
            }

            //  Step 4: Inject context consensus into mediating population.
            //  Assemble the joint candidate from per-party contexts.
            {
                std::vector<double> cand(static_cast<size_t>(total_dim), 0.5);
                for (int dm = 0; dm < num_dms; dm++) {
                    const auto& idx = m_bench->getPartyIndices(dm);
                    const auto& ctx = m_context_vectors[static_cast<size_t>(dm)];
                    for (size_t i = 0; i < idx.size() && i < ctx.size(); i++)
                        cand[static_cast<size_t>(idx[i])] = ctx[i];
                }
                const double cf = m_bench->evaluateJoint(cand, m_env.get());
                auto worst_it = std::min_element(
                    m_mediating_fitness.begin(), m_mediating_fitness.end());
                if (cf > *worst_it) {
                    const size_t wi = static_cast<size_t>(
                        std::distance(m_mediating_fitness.begin(), worst_it));
                    m_mediating_pop[wi] = cand;
                    m_mediating_fitness[wi] = cf;
                }
            }

            // Periodic debug
            if (m_generation % 20 == 0) {
                std::cout << "[gen " << std::setw(4) << m_generation
                    << "] best_joint=" << std::fixed << std::setprecision(4)
                    << getBestConsensusFitness() << " | ctx:";
                for (int dm = 0; dm < num_dms; dm++) {
                    std::cout << " DM" << dm << "=[";
                    for (size_t i = 0;
                        i < m_context_vectors[static_cast<size_t>(dm)].size(); i++) {
                        std::cout << std::setprecision(4)
                            << m_context_vectors[static_cast<size_t>(dm)][i];
                        if (i + 1 < m_context_vectors[static_cast<size_t>(dm)].size())
                            std::cout << ",";
                    }
                    std::cout << "]";
                }
                std::cout << std::endl;
            }
        }

        double getBestConsensusFitness() const {
            if (m_mediating_fitness.empty())
                return -std::numeric_limits<double>::max();
            return *std::max_element(
                m_mediating_fitness.begin(), m_mediating_fitness.end());
        }

        std::vector<double> getBestConsensusSolution() const {
            if (m_mediating_pop.empty()) return {};
            const auto it = std::max_element(
                m_mediating_fitness.begin(), m_mediating_fitness.end());
            return m_mediating_pop[static_cast<size_t>(
                std::distance(m_mediating_fitness.begin(), it))];
        }

        void saveSolutionsAtGeneration(int gen_num) const {
            const std::string fn =
                m_output_directory + "/gen_" + std::to_string(gen_num) + "_solutions.txt";
            std::ofstream out(fn);
            if (!out.is_open()) return;
            out << std::fixed << std::setprecision(6)
                << "# gen=" << gen_num << "  fmt: x1 x2 ... xd joint_fitness\n";
            for (size_t i = 0; i < m_mediating_pop.size(); i++) {
                for (const double v : m_mediating_pop[i]) out << v << " ";
                out << m_mediating_fitness[i] << "\n";
            }
        }

        const std::vector<std::vector<double>>& getMediatingPopulation() const {
            return m_mediating_pop;
        }

    private:
        void updateContextFromMediating() {
            const int num_dms = m_bench->getNumDMs();
            m_context_vectors.resize(static_cast<size_t>(num_dms));

            if (m_mediating_pop.empty()) {
                for (int dm = 0; dm < num_dms; dm++)
                    m_context_vectors[static_cast<size_t>(dm)].assign(
                        static_cast<size_t>(m_bench->getPartyDim(dm)), 0.5);
                return;
            }

            const size_t best_i = static_cast<size_t>(std::distance(
                m_mediating_fitness.begin(),
                std::max_element(m_mediating_fitness.begin(),
                    m_mediating_fitness.end())));
            const auto& best_sol = m_mediating_pop[best_i];

            for (int dm = 0; dm < num_dms; dm++) {
                const auto& idx = m_bench->getPartyIndices(dm);
                const int   dm_dim = m_bench->getPartyDim(dm);
                m_context_vectors[static_cast<size_t>(dm)].resize(
                    static_cast<size_t>(dm_dim));
                for (int d = 0; d < dm_dim; d++)
                    m_context_vectors[static_cast<size_t>(dm)][static_cast<size_t>(d)] =
                    best_sol[static_cast<size_t>(idx[d])];
            }
        }
    };

    //  CEC2015 Metrics
    struct CEC2015Metrics {
        double SR = 0.0, ANOF = 0.0, SP = 0.0, MPR = 0.0;
    };

    CEC2015Metrics evaluateCEC2015Metrics(
        const std::vector<std::vector<double>>& candidates,
        FreePeaks* problem,
        Environment* env,
        double eps_x = 0.1,
        double eps_f = 1.0)
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
            ? 200.0 * problem->numberVariables() * n_opt
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

    //  createComplex2DProblem
    //
    //  Peak centers are specified in the FUNCTION DOMAIN
    //  [-100, 100]^2, NOT in [0,1] space.
    //
    //  The FunctionBase::evaluate() remaps a point from the subspace [a,b]^d
    //  to [-100, 100]^d before computing distances to peak centers.  If you
    //  set a center to {0.3, 0.7} it sits at (0.3, 0.7) in [-100,100] space —
    //  essentially AT the origin — so any point in the domain is within ~1 unit
    //  of the peak, making all evaluations return near-maximum fitness and
    //  removing any landscape structure.
    //
    //  With centers in [-100, 100] space, e.g. {0, 0}, the subspace midpoint
    //  maps exactly there, placing the peak at the center of each subspace.

    void createComplex2DProblem(const std::string& problem_name = "complex_2d") {
        using namespace ofec;
        using namespace free_peaks;

        const std::string problem_dir = "multiparty_multimodal/";
        const std::string vis_dir = "visualization/";

        std::filesystem::create_directories(
            g_working_directory + "/instance/problem/continuous/free_peaks/" + problem_dir);
        std::filesystem::create_directories(
            g_working_directory + "/instance/problem/continuous/free_peaks/" + vis_dir);
        std::filesystem::create_directories(
            g_working_directory + "/instance/problem/continuous/free_peaks/subproblem/" + vis_dir);
        std::filesystem::create_directories(
            g_working_directory + "/instance/problem/continuous/free_peaks/subproblem/function/one_peak/" + vis_dir);

        const int numDim = 2, numObj = 1, numCon = 0;
        std::shared_ptr<ofec::Random> rnd(new Random(0.42));

        FreePeaks::registerFP();
        const std::string freepeakName = "free_peaks";
        std::shared_ptr<Environment> env(generateEnvironmentByFactory(freepeakName));
        env->recordInputParameters();
        env->initialize();
        env->setProblem(generateProblemByFactory(freepeakName));
        auto* freepeak = dynamic_cast<FreePeaks*>(env->problem());

        ParameterMap fp;
        fp["generation_type"] = std::string("assigned");
        fp["dataFile1"] = problem_dir + problem_name + ".txt";
        freepeak->inputParameters().input(fp);
        freepeak->initialize(env.get());
        freepeak->setRandom(rnd);
        freepeak->setSizes(numDim, numObj, numCon);

        freepeak->setKDtree({
            {"root", { {"dm0_subspace", 0.5}, {"dm1_subspace", 0.5} }}
            });

        //  Peaks in FUNCTION DOMAIN [-100, 100]^2.
        //
        //  The subspace midpoint (0.25, 0.5) for dm0 maps to (0, 0).
        //  The subspace midpoint (0.75, 0.5) for dm1 also maps to (0, 0)
        //  (each subspace independently maps its range to [-100, 100]).
        //
        //  A peak at function-space (0, 0) is therefore accessible from
        //  BOTH dm0 and dm1 (each maps its midpoint to it).
        //
        //  Peaks at off-center positions create varied landscape structure
        struct PeakSpec {
            double fx, fy; // function-domain coordinates in [-100, 100]
            double height;
            std::string shape; // "s1" linear, "s2" exponential, etc.
        };
        // Global optimum at (0, 0): symmetric, accessible from all subspaces.
        // Local optima at various positions with decreasing heights.
        const std::vector<PeakSpec> peaks = {
            {   0,    0, 90.0, "s1" },  // global optimum
            {  40,   30, 70.0, "s2" },  // local: upper-right
            { -40,  -30, 60.0, "s1" },  // local: lower-left
            {  50,  -40, 50.0, "s3" },  // local: lower-right
            { -50,   40, 45.0, "s2" },  // local: upper-left
            {  20,  -60, 35.0, "s1" },  // local: far lower-right
            { -20,   60, 30.0, "s3" },  // local: far upper-left
        };

        // Party-specific transform parameters (for landscape differentiation)
        struct PartySpec {
            double rotation; double bias; double condition;
        };
        const std::vector<PartySpec> party_specs = {
            { std::numbers::pi / 8, 0.03,  50.0 },  // DM0: mild rotation, low conditioning
            { std::numbers::pi / 6, 0.04, 100.0 },  // DM1: steeper rotation, higher conditioning
        };

        const std::vector<std::string> subspaces = { "dm0_subspace", "dm1_subspace" };

        for (size_t s = 0; s < subspaces.size(); s++) {
            const std::string& sname = subspaces[s];
            const auto& ps = party_specs[s];

            ParameterMap subpro_param;
            subpro_param["subspace"] = sname;
            subpro_param["generation_type"] = std::string("assigned");
            subpro_param["dataFile1"] = vis_dir + sname + ".txt";
            auto subpro(Subproblem::create());
            subpro->initialize(subpro_param, freepeak);

            // Euclidean distance
            {
                auto dis(FactoryFP<DistanceBase>::produce("Euclidean"));
                ParameterMap dp;
                dis->initialize(freepeak, sname, dp);
                subpro->setDistance(dis);
            }

            // One-peak function
            {
                ParameterMap fun_p;
                fun_p["generation_type"] = std::string("assigned");
                fun_p["dataFile1"] = vis_dir + sname + "_onepeak.txt";
                auto func(FactoryFP<FunctionBase>::produce("one_peak"));
                func->initialize(freepeak, sname, fun_p);
                auto* opf = dynamic_cast<OnePeakFunction*>(func);

                for (const auto& pk : peaks) {
                    auto onepeak(FactoryFP<OnePeakBase>::produce(pk.shape));
                    ParameterMap op;
                    op["center_type"] = std::string("assigned");
                    op["height"] = pk.height;

                    const double bias_x =
                        ps.bias * std::sin(static_cast<double>(s) * 2.1) * 200.0;
                    const double bias_y =
                        ps.bias * std::cos(static_cast<double>(s) * 1.9) * 200.0;

                    std::vector<Real> center = {
                        static_cast<Real>(pk.fx + bias_x),
                        static_cast<Real>(pk.fy + bias_y)
                    };
                    // Clamp to valid function domain
                    for (auto& c : center) {
                        if (c < -100.f) c = -100.f;
                        if (c > 100.f) c = 100.f;
                    }
                    op["center_postion"] = center;
                    onepeak->initialize(freepeak, sname, op);
                    opf->addOnePeaks(onepeak);
                }
                subpro->setFunction(func);
            }

            // Party-bias X transform
            {
                auto t(FactoryFP<X_TransformBase>::produce("MapXPartyBias"));
                ParameterMap tp;
                tp["party_id"] = static_cast<int>(s);
                tp["magnitude"] = ps.bias;
                tp["rotation_angle"] = ps.rotation;
                tp["condition_number"] = ps.condition;
                t->initialize(freepeak, sname, tp);
                subpro->addVariableTransform(t);
            }

            // Asymmetry X transform
            {
                auto t(FactoryFP<X_TransformBase>::produce("MapXAssymetrix"));
                ParameterMap tp;
                tp["alpha"] = Real(0.3 + static_cast<double>(s) * 0.2);
                t->initialize(freepeak, sname, tp);
                subpro->addVariableTransform(t);
            }

            // Objective map
            {
                auto ot(FactoryFP<Y_TransformBase>::produce("map_objective"));
                ParameterMap otp;
                ot->initialize(freepeak, sname, otp);
                subpro->addObjectiveTransform(ot);
            }

            freepeak->setSubproblem(sname, subpro);
        }

        freepeak->bindData();

        // Persist to file
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
            << " peaks in [-100,100]^2 function space, 2 parties." << std::endl;
    }


    //  outputComplexLandscapes
    void outputComplexLandscapes(Environment* env, const std::string& output_dir) {
        std::filesystem::create_directories(output_dir);
        const int res = 200;
        std::cout << ">>> Sampling landscape (" << res << "x" << res << ")..." << std::flush;

        // Conditional slices
        for (double x2_fixed : {0.25, 0.50, 0.75}) {
            std::string tag = std::to_string(x2_fixed).substr(0, 4);
            std::ofstream f(output_dir + "/slice_x2_" + tag + ".txt");
            f << std::fixed << std::setprecision(6);
            for (int i = 0; i < res; i++) {
                const double x1 = static_cast<double>(i) / (res - 1);
                auto sol = env->problem()->createSolution();
                auto& v = dynamic_cast<Solution<VariableVector<Real>>&>(*sol)
                    .variable().vector();
                v[0] = static_cast<Real>(x1);
                v[1] = static_cast<Real>(x2_fixed);
                sol->evaluate(env, false);
                f << x1 << " " << sol->objective(0) << "\n";
            }
        }

        // Full grid
        std::ofstream gf(output_dir + "/global_grid.txt");
        gf << std::fixed << std::setprecision(5);
        for (int i = 0; i < res; i++) {
            for (int j = 0; j < res; j++) {
                const double x1 = static_cast<double>(i) / (res - 1);
                const double x2 = static_cast<double>(j) / (res - 1);
                auto sol = env->problem()->createSolution();
                auto& v = dynamic_cast<Solution<VariableVector<Real>>&>(*sol)
                    .variable().vector();
                v[0] = static_cast<Real>(x1); v[1] = static_cast<Real>(x2);
                sol->evaluate(env, false);
                gf << x1 << " " << x2 << " " << sol->objective(0) << "\n";
            }
        }
        std::cout << " done." << std::endl;
    }

    //  runMPMCoEAExperiment
    void runMPMCoEAExperiment() {
        ofec::registerInstance();

        ofec::g_working_directory =
            R"(E:\HITSZ\Research\Multimodal_Multiparty_Optimization\ThesisProject\Data\ofec_data_new)";

        const std::string vis_base =
            R"(E:\HITSZ\Research\Multimodal_Multiparty_Optimization\ThesisProject\Visualization)";

        std::cout << "\n>>> [MPM-CoEA] Creating benchmark..." << std::endl;
        createComplex2DProblem("complex_2d");

        // Load
        const std::string fpname = "free_peaks";
        auto* env_raw = generateEnvironmentByFactory(fpname);
        env_raw->recordInputParameters();
        env_raw->initialize();
        env_raw->setProblem(generateProblemByFactory(fpname));

        ParameterMap params;
        params["generation_type"] = std::string("read_file");
        params["dataFile1"] = std::string("multiparty_multimodal/complex_2d.txt");
        env_raw->problem()->inputParameters().input(params);
        env_raw->problem()->recordInputParameters();
        env_raw->initializeProblem(0.5);

        auto* fp = CAST_FPs(env_raw->problem());
        if (!fp) throw std::runtime_error("Problem is not FreePeaks.");

        std::cout << ">>> Loaded: dim=" << fp->numberVariables()
            << " | parties=2 | peaks=7" << std::endl;

        outputComplexLandscapes(env_raw, vis_base + "/landscape_before");

        auto bench = std::make_shared<MPMMO_Benchmark>(fp, 2, env_raw);
        auto env_sp = std::shared_ptr<Environment>(env_raw, [](Environment*) {});

        MPM_CoEA solver(bench, env_sp,
            /*pop_size*/   30,
            /*med_pop*/   200,
            /*K*/           5,
            vis_base + "/solutions");

        solver.initialize();
        solver.saveSolutionsAtGeneration(0);
        std::cout << "Gen   0 | " << std::fixed << std::setprecision(4)
            << solver.getBestConsensusFitness() << std::endl;

        const int MAX_GEN = 200;
        for (int gen = 1; gen <= MAX_GEN; gen++) {
            solver.run_one_generation();
            if (gen %  5 == 0 || gen == MAX_GEN) {
                std::cout << "Gen " << std::setw(4) << gen << " | "
                    << std::fixed << std::setprecision(4)
                    << solver.getBestConsensusFitness() << std::endl;
            }
            if (gen % 50 == 0 || gen == MAX_GEN)
                solver.saveSolutionsAtGeneration(gen);
        }

        // CEC2015 metrics
        auto metrics = evaluateCEC2015Metrics(
            solver.getMediatingPopulation(), fp, env_raw);
        std::cout << "\n>>> CEC2015 Metrics\n"
            << "    SR   : " << metrics.SR << "%\n"
            << "    ANOF : " << metrics.ANOF << "\n"
            << "    MPR  : " << metrics.MPR << "\n"
            << "    SP   : " << metrics.SP << std::endl;

        outputComplexLandscapes(env_raw, vis_base + "/landscape_after");
        std::cout << "\n>>> [MPM-CoEA] Complete." << std::endl;
    }

} // namespace ofec
#endif // OFEC_MPMCOEA_SOLVER_HPP