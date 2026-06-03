#ifndef OFEC_MPMCOEA_SOLVER_HPP
#define OFEC_MPMCOEA_SOLVER_HPP

#include "custom_method.hpp"
#include "../instance/problem/continuous/free_peaks/free_peaks.h"
#include "../instance/problem/continuous/free_peaks/subproblem/transform/transform_x/transform_x_base.h"
#include "../instance/problem/continuous/free_peaks/subproblem/transform/transform_y/transform_y_base.h"
#include "../instance/problem/continuous/free_peaks/factory.h"
#include "../instance/problem/continuous/free_peaks/subproblem/subproblem.h"
#include "../instance/problem/continuous/free_peaks/subproblem/function/one_peak_function.h"
#include "../core/environment/environment.h"

#include <queue>
#include <numeric>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <memory>
#include <random>
#include <filesystem>
#include <string>
#include <limits>
#include <sstream>
#include <functional>
#include <set>
#include <map>

namespace ofec {

    //  MPMMO_Benchmark
    class MPMMO_Benchmark {
    private:
        FreePeaks* m_problem;
        Environment* m_env;
        int          m_num_dms;
        std::vector<std::vector<int>> m_party_indices;

    public:
        MPMMO_Benchmark(FreePeaks* prob, int num_dms, Environment* env)
            : m_problem(prob), m_env(env), m_num_dms(num_dms) {
            decomposeVariables();
        }

        void decomposeVariables() {
            const int D = static_cast<int>(m_problem->numberVariables());
            m_party_indices.resize(static_cast<size_t>(m_num_dms));
            for (int v = 0; v < D; v++)
                m_party_indices[static_cast<size_t>(v % m_num_dms)].push_back(v);
            std::cout << "[INFO] Variable decomposition (dim=" << D
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

        double evaluateJoint(const std::vector<double>& x, Environment* env) const {
            std::unique_ptr<SolutionBase> sol(m_problem->createSolution());
            auto& sv = dynamic_cast<Solution<VariableVector<Real>>&>(*sol).variable().vector();
            for (size_t i = 0; i < sv.size() && i < x.size(); i++) sv[i] = static_cast<Real>(x[i]);
            sol->evaluate(env, false);
            double r = static_cast<double>(sol->objective(0));
            return std::isfinite(r) ? r : -1e30;
        }

        double evaluateDM(int dm, const std::vector<double>& my_vars,
            const std::vector<std::vector<double>>& ctx, Environment* env) const {
            std::vector<double> full(static_cast<size_t>(m_problem->numberVariables()), 0.5);
            for (size_t i = 0; i < m_party_indices[dm].size() && i < my_vars.size(); i++)
                full[static_cast<size_t>(m_party_indices[dm][i])] = my_vars[i];
            for (int other = 0; other < m_num_dms; other++) {
                if (other == dm) continue;
                for (size_t i = 0; i < m_party_indices[other].size() && i < ctx[other].size(); i++)
                    full[static_cast<size_t>(m_party_indices[other][i])] = ctx[other][i];
            }
            return evaluateJoint(full, env);
        }

        int getNumDMs() const { return m_num_dms; }
        int getNumVariables() const { return static_cast<int>(m_problem->numberVariables()); }
        int getPartyDim(int p) const { return static_cast<int>(m_party_indices[static_cast<size_t>(p)].size()); }
        const  std::vector<int>& getPartyIndices(int p) const { return m_party_indices[static_cast<size_t>(p)]; }
        FreePeaks* getProblem() const { return m_problem; }
        Environment* getEnv() const { return m_env; }
    };

    //  Nearest-Better Clustering
    std::vector<std::vector<int>> nearestBetterClustering(
        const std::vector<std::vector<double>>& pop,
        const std::vector<double>& fit,
        double phi = 2.0)
    {
        const int n = static_cast<int>(pop.size());
        if (n == 0) return {};
        const int dim = static_cast<int>(pop[0].size());

        struct Edge { int from, to; double w; };
        std::vector<Edge> edges;
        edges.reserve(static_cast<size_t>(n));
        for (int i = 0; i < n; i++) {
            double best_d = std::numeric_limits<double>::max();
            int    best_j = -1;
            for (int j = 0; j < n; j++) {
                if (fit[j] <= fit[i] + 1e-12) continue;
                double d2 = 0.0;
                for (int k = 0; k < dim; k++) { double dv = pop[i][k] - pop[j][k]; d2 += dv * dv; }
                if (std::sqrt(d2) < best_d) { best_d = std::sqrt(d2); best_j = j; }
            }
            if (best_j >= 0) edges.push_back({ i, best_j, best_d });
        }

        if (edges.empty()) {
            std::vector<std::vector<int>> cl; cl.reserve(static_cast<size_t>(n));
            for (int i = 0; i < n; i++) cl.push_back({ i });
            return cl;
        }

        double sum_w = 0.0;
        for (const auto& e : edges) sum_w += e.w;
        const double thr = phi * sum_w / static_cast<double>(edges.size());

        std::vector<std::vector<int>> adj(static_cast<size_t>(n));
        for (const auto& e : edges)
            if (e.w <= thr) { adj[e.from].push_back(e.to); adj[e.to].push_back(e.from); }

        std::vector<bool> visited(static_cast<size_t>(n), false);
        std::vector<std::vector<int>> clusters;
        for (int i = 0; i < n; i++) {
            if (visited[i]) continue;
            std::vector<int> comp;
            std::queue<int> q; q.push(i); visited[i] = true;
            while (!q.empty()) {
                int u = q.front(); q.pop(); comp.push_back(u);
                for (int v : adj[u]) if (!visited[v]) { visited[v] = true; q.push(v); }
            }
            clusters.push_back(std::move(comp));
        }
        return clusters;
    }

    //  CEC-2015 Metrics
    struct CEC2015Metrics { double SR = 0, ANOF = 0, SP = 0, MPR = 0; };

    CEC2015Metrics evaluateCEC2015Metrics(
        const std::vector<std::vector<double>>& candidates,
        FreePeaks* problem, Environment* env,
        double eps_x = 0.05, double eps_f = 2.0)
    {
        CEC2015Metrics m;
        const auto* opt = problem->optima();
        if (!opt || opt->numberSolutions() == 0) return m;
        const int nopt = static_cast<int>(opt->numberSolutions());
        std::vector<bool> found(static_cast<size_t>(nopt), false);
        int total_found = 0;
        for (const auto& cand : candidates) {
            std::unique_ptr<SolutionBase> sol(problem->createSolution());
            auto& sv = dynamic_cast<Solution<VariableVector<Real>>&>(*sol).variable().vector();
            for (size_t d = 0; d < cand.size() && d < sv.size(); d++) sv[d] = static_cast<Real>(cand[d]);
            sol->evaluate(env, false);
            const double cf = sol->objective(0);
            if (!std::isfinite(cf)) continue;
            for (int i = 0; i < nopt; i++) {
                if (found[i]) continue;
                const auto& ov = opt->solution(i).variable().vector();
                double d2 = 0.0;
                for (size_t d = 0; d < cand.size() && d < ov.size(); d++) {
                    double dv = cand[d] - static_cast<double>(ov[d]); d2 += dv * dv;
                }
                if (std::sqrt(d2) < eps_x && std::abs(cf - opt->solution(i).objective(0)) < eps_f)
                {
                    found[i] = true; ++total_found;
                }
            }
        }
        m.ANOF = static_cast<double>(total_found) / nopt;
        m.SR = (total_found == nopt) ? 100.0 : 0.0;
        m.SP = (m.SR > 0.0) ? static_cast<double>(candidates.size()) * nopt
            : std::numeric_limits<double>::infinity();
        double gf = -std::numeric_limits<double>::max();
        for (int i = 0; i < nopt; i++) gf = std::max(gf, opt->solution(i).objective(0));
        double sn = 0, sd = 0;
        for (int i = 0; i < nopt; i++) {
            double w = (gf - opt->solution(i).objective(0)) + 1.0;
            if (found[i]) sn += w; sd += w;
        }
        m.MPR = (sd > 0) ? sn / sd : 0.0;
        return m;
    }

    //  MPM-CoEA
    class MPM_CoEA {
    private:
        std::shared_ptr<MPMMO_Benchmark> m_bench;
        std::shared_ptr<Environment>     m_env;

        // Hyperparameters
        int    m_K;
        int    m_niche_pop;
        int    m_med_pop_size;
        int    m_party_gens;
        int    m_coevo_gens;
        int    m_refine_gens;
        double m_sigma;
        double m_de_F;
        double m_de_CR;
        double m_sharing_radius;
        double m_sharing_alpha;

        // Stagnation
        static constexpr int    STAGNATION_LIMIT = 60;
        static constexpr int    MIN_POP = 15;
        static constexpr double DEAD = -1e30;

        // Per-party data
        std::vector<std::vector<std::vector<double>>> m_party_optima;
        std::vector<std::vector<double>>              m_party_opt_fit;
        std::vector<std::vector<std::vector<std::vector<double>>>> m_competing_pops;
        std::vector<std::vector<std::vector<double>>>              m_competing_fit;

        // Mediating population
        std::vector<std::vector<double>> m_mediating_pop;
        std::vector<double>              m_mediating_fit;

        // DM-oriented elites
        std::vector<std::vector<std::vector<double>>> m_dm_elite;
        std::vector<std::vector<double>>              m_dm_elite_fit;

        // Context archive
        std::vector<std::vector<std::vector<double>>> m_context_archive;
        std::vector<std::vector<double>>              m_context_fit;

        //  Optimal Solution Archive (OSA)
        std::vector<std::vector<double>> m_osa;
        std::vector<double>              m_osa_fit;
        double                           m_osa_radius;
        static constexpr int             OSA_MAX = 300;

        // Diversity-aware OSA insertion
        void osaInsert(const std::vector<double>& sol, double f) {
            if (f <= DEAD + 1.0 || !std::isfinite(f)) return;

            // Check distance to all existing OSA entries
            double min_d = std::numeric_limits<double>::max();
            size_t closest_i = 0;
            for (size_t i = 0; i < m_osa.size(); ++i) {
                double d2 = 0.0;
                for (size_t j = 0; j < sol.size(); ++j) {
                    double dv = sol[j] - m_osa[i][j];
                    d2 += dv * dv;
                }
                double d = std::sqrt(d2);
                if (d < min_d) { min_d = d; closest_i = i; }
            }

            if (static_cast<int>(m_osa.size()) < OSA_MAX) {
                if (min_d >= m_osa_radius) {
                    m_osa.push_back(sol); m_osa_fit.push_back(f);
                }
                else if (f > m_osa_fit[closest_i]) {
                    // Same niche but better: replace
                    m_osa[closest_i] = sol; m_osa_fit[closest_i] = f;
                }
            }
            else {
                if (min_d > m_osa_radius) {
                    // New niche discovered: evict weakest incumbent regardless of fitness
                    auto wt = std::min_element(m_osa_fit.begin(), m_osa_fit.end());
                    size_t wi = static_cast<size_t>(std::distance(m_osa_fit.begin(), wt));
                    m_osa[wi] = sol; m_osa_fit[wi] = f;
                }
                else if (f > m_osa_fit[closest_i]) {
                    // Same niche, better fitness
                    m_osa[closest_i] = sol; m_osa_fit[closest_i] = f;
                }
            }
        }

        void adaptOsaRadius() {
            double progress = static_cast<double>(m_generation) / m_coevo_gens;
            m_osa_radius = m_sigma * (1.0 - 0.50 * progress);
            m_osa_radius = std::max(0.03, m_osa_radius);
        }

        // Stagnation
        double m_best_prev = DEAD;
        int    m_stagnation_cnt = 0;
        int    m_restart_count = 0;
        int    m_generation = 0;

        // Output
        std::string m_output_dir;
        std::vector<std::pair<int, double>> m_conv_log;

        // Helpers
        double evalJoint(const std::vector<double>& x) const
        {
            return m_bench->evaluateJoint(x, m_env.get());
        }

        double evalDM(int dm, const std::vector<double>& my,
            const std::vector<std::vector<double>>& ctx) const
        {
            return m_bench->evaluateDM(dm, my, ctx, m_env.get());
        }

        std::vector<std::vector<double>> neutralCtx() const {
            const int nd = m_bench->getNumDMs();
            std::vector<std::vector<double>> ctx(static_cast<size_t>(nd));
            for (int d = 0; d < nd; d++)
                ctx[d].assign(static_cast<size_t>(m_bench->getPartyDim(d)), 0.5);
            return ctx;
        }

        std::vector<std::vector<double>> ctxFromJoint(int dm, const std::vector<double>& x) const {
            const int nd = m_bench->getNumDMs();
            auto ctx = neutralCtx();
            for (int other = 0; other < nd; other++) {
                if (other == dm) continue;
                const auto& idx = m_bench->getPartyIndices(other);
                ctx[other].resize(idx.size());
                for (size_t k = 0; k < idx.size(); k++) ctx[other][k] = x[idx[k]];
            }
            return ctx;
        }

        std::vector<std::vector<double>> lhsSample(int n, int dim, std::mt19937& rng) {
            std::uniform_real_distribution<double> u(0.0, 1.0);
            std::vector<std::vector<double>> pop(n, std::vector<double>(dim));
            for (int d = 0; d < dim; d++) {
                std::vector<double> st(n);
                for (int i = 0; i < n; i++) st[i] = (i + u(rng)) / n;
                std::shuffle(st.begin(), st.end(), rng);
                for (int i = 0; i < n; i++) pop[i][d] = st[i];
            }
            return pop;
        }

        static double clamp(double v) { return std::max(0.0, std::min(1.0, v)); }

        //  Phase 1a: NBC init
        void NBC_init(int dm, const std::vector<std::vector<double>>& ctx) {
            const int dim = m_bench->getPartyDim(dm);
            const int pool = std::max(300, m_K * m_niche_pop * 3);
            std::mt19937 rng(static_cast<unsigned>(dm * 12345u + 1u));
            auto raw_pop = lhsSample(pool, dim, rng);
            std::vector<double> raw_fit(pool);
            for (int i = 0; i < pool; i++) raw_fit[i] = evalDM(dm, raw_pop[i], ctx);

            auto clusters = nearestBetterClustering(raw_pop, raw_fit, 2.0);
            std::sort(clusters.begin(), clusters.end(),
                [&](const std::vector<int>& a, const std::vector<int>& b) {
                    double ba = DEAD, bb = DEAD;
                    for (int i : a) ba = std::max(ba, raw_fit[i]);
                    for (int i : b) bb = std::max(bb, raw_fit[i]);
                    return ba > bb;
                });
            if (static_cast<int>(clusters.size()) > m_K) clusters.resize(m_K);

            m_competing_pops[dm].clear();
            m_competing_fit[dm].clear();
            std::uniform_real_distribution<double> u(0.0, 1.0);

            for (auto& cl : clusters) {
                std::sort(cl.begin(), cl.end(), [&](int a, int b) { return raw_fit[a] > raw_fit[b]; });
                const int take = std::min(static_cast<int>(cl.size()), m_niche_pop);
                std::vector<std::vector<double>> sp;
                std::vector<double> sf;
                for (int k = 0; k < take; k++) { sp.push_back(raw_pop[cl[k]]); sf.push_back(raw_fit[cl[k]]); }
                std::normal_distribution<double> gauss(0.0, 0.05);
                while (static_cast<int>(sp.size()) < MIN_POP) {
                    auto perturbed = sp[0];
                    for (auto& v : perturbed) v = clamp(v + gauss(rng));
                    sp.push_back(perturbed);
                    sf.push_back(evalDM(dm, perturbed, ctx));
                }
                m_competing_pops[dm].push_back(std::move(sp));
                m_competing_fit[dm].push_back(std::move(sf));
            }
            std::cout << "  Party " << dm << ": NBC :: " << m_competing_pops[dm].size() << " sub-pops\n";
        }

        //  Phase 1b: Clearing-DE per party
        void ClearingDE(int dm, const std::vector<std::vector<double>>& ctx,
            unsigned seed = 0) {
            const int dim = m_bench->getPartyDim(dm);
            std::mt19937 rng(seed == 0 ? static_cast<unsigned>(dm * 77777u + 3u) : seed);
            std::uniform_real_distribution<double> u(0.0, 1.0);

            for (size_t sp = 0; sp < m_competing_pops[dm].size(); sp++) {
                auto& pop = m_competing_pops[dm][sp];
                auto& fit = m_competing_fit[dm][sp];
                const int np = static_cast<int>(pop.size());
                if (np < 4) continue;
                std::uniform_int_distribution<int> pick(0, np - 1);

                for (int gen = 0; gen < m_party_gens; gen++) {
                    for (int i = 0; i < np; i++) {
                        if (fit[i] <= DEAD + 1.0) continue;
                        int r1 = i, r2 = i, r3 = i, att = 0;
                        while (r1 == i && ++att < 2000) r1 = pick(rng);
                        att = 0; while ((r2 == i || r2 == r1) && ++att < 2000) r2 = pick(rng);
                        att = 0; while ((r3 == i || r3 == r1 || r3 == r2) && ++att < 2000) r3 = pick(rng);
                        if (r1 == i || r2 == i || r3 == i) continue;

                        std::normal_distribution<double> gF(m_de_F, 0.1);
                        const double F = std::max(0.1, std::min(0.9, gF(rng)));

                        std::vector<double> trial(dim);
                        const int jr = static_cast<int>(u(rng) * dim);
                        for (int d = 0; d < dim; d++) {
                            if (u(rng) < m_de_CR || d == jr)
                                trial[d] = clamp(pop[r1][d] + F * (pop[r2][d] - pop[r3][d]));
                            else
                                trial[d] = pop[i][d];
                        }
                        const double tf = evalDM(dm, trial, ctx);
                        if (tf > fit[i]) { pop[i] = std::move(trial); fit[i] = tf; }
                    }

                    if (gen % 20 == 0) {
                        for (size_t s2 = 0; s2 < m_competing_pops[dm].size(); s2++) {
                            for (int j = 0; j < static_cast<int>(m_competing_pops[dm][s2].size()); j++) {
                                if (m_competing_fit[dm][s2][j] <= DEAD + 1.0) continue;
                                for (int i = 0; i < np; i++) {
                                    if (fit[i] <= DEAD + 1.0 || (s2 == sp && i == j)) continue;
                                    if (s2 == sp && fit[i] > m_competing_fit[dm][s2][j]) continue;
                                    double d2 = 0.0;
                                    for (int d = 0; d < dim; d++) {
                                        double dv = pop[i][d] - m_competing_pops[dm][s2][j][d];
                                        d2 += dv * dv;
                                    }
                                    if (std::sqrt(d2) < m_sigma) {
                                        if (fit[i] <= m_competing_fit[dm][s2][j])
                                            fit[i] = DEAD;
                                    }
                                }
                            }
                        }
                    }
                }

                auto bit = std::max_element(fit.begin(), fit.end());
                size_t bi = static_cast<size_t>(std::distance(fit.begin(), bit));
                if (fit[bi] > DEAD + 1.0) {
                    m_party_optima[dm].push_back(pop[bi]);
                    m_party_opt_fit[dm].push_back(fit[bi]);
                }
            }
        }

        //  Phase 2: Cartesian seed
        void CartesianSeed() {
            const int nd = m_bench->getNumDMs();
            const int tot = m_bench->getNumVariables();

            std::vector<std::vector<double>> candidates;
            std::vector<double>              cfit;

            std::vector<int> cur(static_cast<size_t>(nd), 0);
            std::function<void(int)> dfs = [&](int dm) {
                if (dm == nd) {
                    std::vector<double> full(static_cast<size_t>(tot), 0.5);
                    for (int d = 0; d < nd; d++) {
                        if (m_party_optima[d].empty()) return;
                        const auto& idx = m_bench->getPartyIndices(d);
                        const auto& pos = m_party_optima[d][cur[d]];
                        for (size_t k = 0; k < idx.size() && k < pos.size(); k++)
                            full[idx[k]] = pos[k];
                    }
                    candidates.push_back(full);
                    cfit.push_back(evalJoint(full));
                    return;
                }
                for (int i = 0; i < static_cast<int>(m_party_optima[dm].size()); i++)
                {
                    cur[dm] = i; dfs(dm + 1);
                }
                };
            dfs(0);

            std::vector<size_t> ord(candidates.size()); std::iota(ord.begin(), ord.end(), 0);
            std::sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return cfit[a] > cfit[b]; });

            const double min_sep = 0.03;
            m_mediating_pop.clear(); m_mediating_fit.clear();
            for (size_t idx : ord) {
                if (static_cast<int>(m_mediating_pop.size()) >= m_med_pop_size) break;
                bool close = false;
                for (const auto& ex : m_mediating_pop) {
                    double d2 = 0;
                    for (size_t d = 0; d < candidates[idx].size(); d++) {
                        double dv = candidates[idx][d] - ex[d]; d2 += dv * dv;
                    }
                    if (std::sqrt(d2) < min_sep) { close = true; break; }
                }
                if (!close) {
                    m_mediating_pop.push_back(candidates[idx]);
                    m_mediating_fit.push_back(cfit[idx]);
                    osaInsert(candidates[idx], cfit[idx]);
                }
            }
            std::mt19937 rng(123u);
            std::uniform_real_distribution<double> u(0.0, 1.0);
            while (static_cast<int>(m_mediating_pop.size()) < m_med_pop_size) {
                std::vector<double> rsol(tot);
                for (auto& v : rsol) v = u(rng);
                m_mediating_pop.push_back(rsol);
                m_mediating_fit.push_back(evalJoint(rsol));
            }
            std::cout << "  Mediating pop seeded: " << m_mediating_pop.size()
                << " | best=" << std::fixed << std::setprecision(4)
                << *std::max_element(m_mediating_fit.begin(), m_mediating_fit.end()) << "\n";
        }

        //  Phase 3a: Evolve competing sub-pops
        void evolveCompeting() {
            const int nd = m_bench->getNumDMs();
            std::mt19937 rng(static_cast<unsigned>(m_generation * 123457u));
            std::uniform_real_distribution<double> u(0.0, 1.0);

            for (int dm = 0; dm < nd; dm++) {
                const int dim = m_bench->getPartyDim(dm);
                for (size_t sp = 0; sp < m_competing_pops[dm].size(); sp++) {
                    auto& pop = m_competing_pops[dm][sp];
                    auto& fit = m_competing_fit[dm][sp];
                    const int np = static_cast<int>(pop.size());
                    if (np < 4) continue;

                    std::vector<std::vector<double>> ctx = neutralCtx();
                    for (int other = 0; other < nd; other++) {
                        if (other == dm || m_context_archive[other].empty()) continue;
                        size_t ctx_i = sp % m_context_archive[other].size();
                        ctx[other] = m_context_archive[other][ctx_i];
                    }

                    std::uniform_int_distribution<int> pick(0, np - 1);
                    std::vector<std::vector<double>> new_pop = pop;
                    std::vector<double>              new_fit = fit;

                    for (int i = 0; i < np; i++) {
                        if (fit[i] <= DEAD + 1.0) continue;
                        int r1 = i, r2 = i, r3 = i, att = 0;
                        while (r1 == i && ++att < 200)r1 = pick(rng);
                        att = 0; while ((r2 == i || r2 == r1) && ++att < 200)r2 = pick(rng);
                        att = 0; while ((r3 == i || r3 == r1 || r3 == r2) && ++att < 200)r3 = pick(rng);
                        if (r1 == i || r2 == i || r3 == i) continue;

                        std::vector<double> trial(dim);
                        const int jr = static_cast<int>(u(rng) * dim);
                        for (int d = 0; d < dim; d++) {
                            if (u(rng) < m_de_CR || d == jr)
                                trial[d] = clamp(pop[r1][d] + m_de_F * (pop[r2][d] - pop[r3][d]));
                            else trial[d] = pop[i][d];
                        }
                        double tf = evalDM(dm, trial, ctx);
                        if (tf > fit[i]) { new_pop[i] = trial; new_fit[i] = tf; }
                    }

                    // CMA-ES Gaussian mutation on best
                    {
                        std::vector<double> mean(dim, 0.0);
                        int alive = 0;
                        for (int i = 0; i < np; i++) {
                            if (new_fit[i] <= DEAD + 1.0) continue;
                            for (int d = 0; d < dim; d++) mean[d] += new_pop[i][d];
                            alive++;
                        }
                        if (alive > 1) {
                            for (auto& v : mean) v /= alive;
                            double sigma_cma = 0.0;
                            for (int i = 0; i < np; i++) {
                                if (new_fit[i] <= DEAD + 1.0) continue;
                                for (int d = 0; d < dim; d++) {
                                    double dv = new_pop[i][d] - mean[d]; sigma_cma += dv * dv;
                                }
                            }
                            sigma_cma = std::sqrt(sigma_cma / (alive * dim));
                            sigma_cma = std::max(0.005, sigma_cma * 0.5);

                            auto best_it = std::max_element(new_fit.begin(), new_fit.end());
                            size_t bi = static_cast<size_t>(std::distance(new_fit.begin(), best_it));
                            std::normal_distribution<double> gauss(0.0, sigma_cma);
                            for (int trial_i = 0; trial_i < 5; trial_i++) {
                                std::vector<double> mut = new_pop[bi];
                                for (auto& v : mut) v = clamp(v + gauss(rng));
                                double mf = evalDM(dm, mut, ctx);
                                if (mf > new_fit[bi]) { new_pop[bi] = mut; new_fit[bi] = mf; }
                            }
                        }
                    }
                    pop = std::move(new_pop);
                    fit = std::move(new_fit);
                }
            }
        }

        //  Phase 3b: Update context archives
        void updateContextArchives() {
            const int nd = m_bench->getNumDMs();
            for (int dm = 0; dm < nd; dm++) {
                std::vector<std::vector<double>> cands;
                std::vector<double>              cfits;
                for (size_t sp = 0; sp < m_competing_pops[dm].size(); sp++) {
                    auto& f = m_competing_fit[dm][sp];
                    auto it = std::max_element(f.begin(), f.end());
                    if (*it <= DEAD + 1.0) continue;
                    size_t bi = static_cast<size_t>(std::distance(f.begin(), it));
                    cands.push_back(m_competing_pops[dm][sp][bi]);
                    cfits.push_back(*it);
                }
                const auto& pm_idx = m_bench->getPartyIndices(dm);
                std::vector<size_t> pm_ord(m_mediating_pop.size());
                std::iota(pm_ord.begin(), pm_ord.end(), 0);
                std::sort(pm_ord.begin(), pm_ord.end(),
                    [&](size_t a, size_t b) { return m_mediating_fit[a] > m_mediating_fit[b]; });
                for (int r = 0; r < std::min(m_K, static_cast<int>(pm_ord.size())); r++) {
                    std::vector<double> proj(pm_idx.size());
                    for (size_t k = 0; k < pm_idx.size(); k++)
                        proj[k] = m_mediating_pop[pm_ord[r]][pm_idx[k]];
                    cands.push_back(proj);
                    cfits.push_back(m_mediating_fit[pm_ord[r]]);
                }
                std::vector<size_t> ord(cands.size()); std::iota(ord.begin(), ord.end(), 0);
                std::sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return cfits[a] > cfits[b]; });
                m_context_archive[dm].clear(); m_context_fit[dm].clear();
                const double min_sep = 0.10;
                for (size_t idx : ord) {
                    if (static_cast<int>(m_context_archive[dm].size()) >= m_K * 2) break;
                    bool close = false;
                    for (const auto& ex : m_context_archive[dm]) {
                        double d2 = 0;
                        for (size_t k = 0; k < cands[idx].size() && k < ex.size(); k++) {
                            double dv = cands[idx][k] - ex[k]; d2 += dv * dv;
                        }
                        if (std::sqrt(d2) < min_sep) { close = true; break; }
                    }
                    if (!close) {
                        m_context_archive[dm].push_back(cands[idx]);
                        m_context_fit[dm].push_back(cfits[idx]);
                    }
                }
            }
        }

        //  Phase 3c: Evolve mediating population
        void evolveMediating() {
            const int nd = m_bench->getNumDMs();
            const int tot = m_bench->getNumVariables();
            const int npm = static_cast<int>(m_mediating_pop.size());
            if (npm < 8) return;

            std::mt19937 rng(static_cast<unsigned>(m_generation * 444333u));
            std::uniform_real_distribution<double> u(0.0, 1.0);

            // Inject best from competing sub-pops with MULTIPLE contexts
            for (int dm = 0; dm < nd; dm++) {
                for (size_t sp = 0; sp < m_competing_pops[dm].size(); sp++) {
                    auto& f = m_competing_fit[dm][sp];
                    auto it = std::max_element(f.begin(), f.end());
                    if (*it <= DEAD + 1.0) continue;
                    size_t bi = static_cast<size_t>(std::distance(f.begin(), it));
                    const auto& idx = m_bench->getPartyIndices(dm);

                    int other = (dm == 0) ? 1 : 0;
                    int n_ctx = 1;
                    if (other < nd && !m_context_archive[other].empty())
                        n_ctx = static_cast<int>(m_context_archive[other].size());
                    n_ctx = std::min(n_ctx, 5); // limit cost

                    for (int ci = 0; ci < n_ctx; ++ci) {
                        auto ctx = neutralCtx();
                        if (other < nd && !m_context_archive[other].empty())
                            ctx[other] = m_context_archive[other][ci];

                        std::vector<double> full(tot, 0.5);
                        for (size_t k = 0; k < idx.size(); k++)
                            full[idx[k]] = m_competing_pops[dm][sp][bi][k];
                        for (int o = 0; o < nd; o++) {
                            if (o == dm) continue;
                            const auto& oidx = m_bench->getPartyIndices(o);
                            for (size_t k = 0; k < oidx.size() && k < ctx[o].size(); k++)
                                full[oidx[k]] = ctx[o][k];
                        }
                        double jf = evalJoint(full);
                        osaInsert(full, jf);
                        auto wt = std::min_element(m_mediating_fit.begin(), m_mediating_fit.end());
                        if (jf > *wt) {
                            size_t wi = static_cast<size_t>(std::distance(m_mediating_fit.begin(), wt));
                            m_mediating_pop[wi] = full; m_mediating_fit[wi] = jf;
                        }
                    }
                }
            }

            // Species-based fitness-sharing DE
            std::uniform_int_distribution<int> pick(0, npm - 1);
            for (int i = 0; i < npm; i++) {
                int r1 = i, r2 = i, r3 = i, att = 0;
                while (r1 == i && ++att < 500)r1 = pick(rng);
                att = 0; while ((r2 == i || r2 == r1) && ++att < 500)r2 = pick(rng);
                att = 0; while ((r3 == i || r3 == r1 || r3 == r2) && ++att < 500)r3 = pick(rng);
                if (r1 == i || r2 == i || r3 == i)continue;

                std::vector<double> trial(tot);
                const int jr = static_cast<int>(u(rng) * tot);
                for (int d = 0; d < tot; d++) {
                    if (u(rng) < m_de_CR || d == jr)
                        trial[d] = clamp(m_mediating_pop[r1][d] + m_de_F * (m_mediating_pop[r2][d] - m_mediating_pop[r3][d]));
                    else trial[d] = m_mediating_pop[i][d];
                }
                double raw_tf = evalJoint(trial);

                // Fitness sharing
                double shared_trial = raw_tf, shared_cur = m_mediating_fit[i];
                {
                    auto sharingPenalty = [&](const std::vector<double>& sol, double raw) {
                        double nc = 1.0;
                        for (size_t j = 0; j < m_mediating_pop.size(); j++) {
                            if (m_mediating_fit[j] < 20.0) continue;
                            double d2 = 0;
                            for (int d = 0; d < tot; d++) { double dv = sol[d] - m_mediating_pop[j][d]; d2 += dv * dv; }
                            double dist = std::sqrt(d2);
                            if (dist < m_sharing_radius)
                                nc += std::max(0.0, 1.0 - std::pow(dist / m_sharing_radius, m_sharing_alpha));
                        }
                        return raw / nc;
                        };
                    shared_trial = sharingPenalty(trial, raw_tf);
                    shared_cur = sharingPenalty(m_mediating_pop[i], m_mediating_fit[i]);
                }
                if (shared_trial >= shared_cur) {
                    m_mediating_pop[i] = std::move(trial);
                    m_mediating_fit[i] = raw_tf;
                    osaInsert(m_mediating_pop[i], m_mediating_fit[i]);
                }
            }
        }

        //  Phase 3d: DM-oriented preservation
        void applyDMPreservation() {
            const int nd = m_bench->getNumDMs();
            const int npm = static_cast<int>(m_mediating_pop.size());
            for (int dm = 0; dm < nd; dm++) {
                std::vector<std::pair<double, int>> dm_scores;
                for (int i = 0; i < npm; i++)
                    dm_scores.push_back({ m_mediating_fit[i], i });
                std::sort(dm_scores.begin(), dm_scores.end(),
                    [](const std::pair<double, int>& a, const std::pair<double, int>& b) { return a.first > b.first; });

                m_dm_elite[dm].clear(); m_dm_elite_fit[dm].clear();
                const double min_sep = m_sigma;
                for (auto& [sc, idx] : dm_scores) {
                    if (static_cast<int>(m_dm_elite[dm].size()) >= m_K) break;
                    bool close = false;
                    for (const auto& ex : m_dm_elite[dm]) {
                        double d2 = 0;
                        for (int d = 0; d < m_bench->getNumVariables(); d++) {
                            double dv = m_mediating_pop[idx][d] - ex[d]; d2 += dv * dv;
                        }
                        if (std::sqrt(d2) < min_sep * 2.0) { close = true; break; }
                    }
                    if (!close) {
                        m_dm_elite[dm].push_back(m_mediating_pop[idx]);
                        m_dm_elite_fit[dm].push_back(sc);
                        osaInsert(m_mediating_pop[idx], sc);
                    }
                }

                for (const auto& elite : m_dm_elite[dm]) {
                    bool present = false;
                    for (const auto& mpi : m_mediating_pop) {
                        double d2 = 0;
                        for (int d = 0; d < m_bench->getNumVariables(); d++) { double dv = elite[d] - mpi[d]; d2 += dv * dv; }
                        if (std::sqrt(d2) < 1e-6) { present = true; break; }
                    }
                    if (!present) {
                        auto wt = std::min_element(m_mediating_fit.begin(), m_mediating_fit.end());
                        size_t wi = static_cast<size_t>(std::distance(m_mediating_fit.begin(), wt));
                        double ef = evalJoint(elite);
                        if (ef > *wt) { m_mediating_pop[wi] = elite; m_mediating_fit[wi] = ef; }
                    }
                }
            }
        }

        //  Phase 3e: Stagnation + adaptive restart
        bool checkAndRestart() {
            double cur = *std::max_element(m_mediating_fit.begin(), m_mediating_fit.end());
            m_conv_log.push_back({ m_generation, cur });
            if (cur > m_best_prev + 1e-6) {
                m_best_prev = cur; m_stagnation_cnt = 0;
            }
            else { ++m_stagnation_cnt; }

            if (m_stagnation_cnt >= STAGNATION_LIMIT) {
                std::cout << "[RESTART gen=" << m_generation
                    << " best=" << std::fixed << std::setprecision(4) << cur << "]\n";

                // Avoid infinite restarts once global optimum is effectively found
                if (m_restart_count >= 2 && m_best_prev >= 89.0) {
                    m_stagnation_cnt = 0;
                    const int tot = m_bench->getNumVariables();
                    const int n_replace = static_cast<int>(m_mediating_pop.size()) / 10;
                    std::mt19937 rng(static_cast<unsigned>(m_generation * 55555u));
                    std::uniform_real_distribution<double> u(0.0, 1.0);
                    for (int r = 0; r < n_replace; r++) {
                        std::vector<double> rsol(tot);
                        for (auto& v : rsol) v = u(rng);
                        double rf = evalJoint(rsol);
                        auto wt = std::min_element(m_mediating_fit.begin(), m_mediating_fit.end());
                        size_t wi = static_cast<size_t>(std::distance(m_mediating_fit.begin(), wt));
                        if (rf > *wt) { m_mediating_pop[wi] = rsol; m_mediating_fit[wi] = rf; }
                    }
                    return false;
                }

                ++m_restart_count; m_stagnation_cnt = 0;

                const int nd = m_bench->getNumDMs();
                std::mt19937 rng(static_cast<unsigned>(m_generation * 55555u + m_restart_count * 9999u));
                std::uniform_real_distribution<double> u(0.0, 1.0);

                for (int dm = 0; dm < nd; dm++) {
                    auto ctx = neutralCtx();
                    for (size_t sp = 0; sp < m_competing_pops[dm].size(); sp++) {
                        auto& pop = m_competing_pops[dm][sp];
                        const int dim = m_bench->getPartyDim(dm);
                        auto& f = m_competing_fit[dm][sp];
                        auto best_it = std::max_element(f.begin(), f.end());
                        size_t best_sp = static_cast<size_t>(std::distance(f.begin(), best_it));
                        for (size_t k = 0; k < pop.size(); k++) {
                            if (k == best_sp) continue;
                            for (auto& v : pop[k]) v = u(rng);
                            m_competing_fit[dm][sp][k] = evalDM(dm, pop[k], ctx);
                        }
                    }
                }

                for (size_t i = 0; i < m_osa.size(); i++) {
                    auto wt = std::min_element(m_mediating_fit.begin(), m_mediating_fit.end());
                    if (m_osa_fit[i] > *wt) {
                        size_t wi = static_cast<size_t>(std::distance(m_mediating_fit.begin(), wt));
                        m_mediating_pop[wi] = m_osa[i];
                        m_mediating_fit[wi] = m_osa_fit[i];
                    }
                }

                const int tot = m_bench->getNumVariables();
                const int n_replace = static_cast<int>(m_mediating_pop.size()) / 5;
                for (int r = 0; r < n_replace; r++) {
                    std::vector<double> rsol(tot);
                    for (auto& v : rsol) v = u(rng);
                    double rf = evalJoint(rsol);
                    auto wt = std::min_element(m_mediating_fit.begin(), m_mediating_fit.end());
                    size_t wi = static_cast<size_t>(std::distance(m_mediating_fit.begin(), wt));
                    if (rf > *wt) { m_mediating_pop[wi] = rsol; m_mediating_fit[wi] = rf; }
                }
                return true;
            }
            return false;
        }

        //  Phase 4: Context-aware refinement
        void ContextRefine(std::vector<double>& cand, int rounds = 3) {
            const int nd = m_bench->getNumDMs();
            std::mt19937 rng(static_cast<unsigned>(m_generation * 777u));
            std::normal_distribution<double> g(0.0, 0.03);
            std::uniform_real_distribution<double> u(0.0, 1.0);
            for (int round = 0; round < rounds; round++) {
                for (int dm = 0; dm < nd; dm++) {
                    const int ddim = m_bench->getPartyDim(dm);
                    const auto& idx = m_bench->getPartyIndices(dm);
                    auto ctx = ctxFromJoint(dm, cand);
                    std::vector<double> my(ddim);
                    for (int k = 0; k < ddim; k++) my[k] = cand[idx[k]];

                    const int mn = 20, mg = 50;
                    std::vector<std::vector<double>> mp(mn, my);
                    std::vector<double> mf(mn);
                    mf[0] = evalDM(dm, my, ctx);
                    for (int i = 1; i < mn; i++) {
                        for (auto& v : mp[i]) v = clamp(v + g(rng));
                        mf[i] = evalDM(dm, mp[i], ctx);
                    }
                    std::uniform_int_distribution<int> pick(0, mn - 1);
                    for (int gg = 0; gg < mg; gg++) for (int i = 0; i < mn; i++) {
                        int r1 = i, r2 = i, r3 = i, att = 0;
                        while (r1 == i && ++att < 100)r1 = pick(rng);
                        att = 0; while ((r2 == i || r2 == r1) && ++att < 100)r2 = pick(rng);
                        att = 0; while ((r3 == i || r3 == r1 || r3 == r2) && ++att < 100)r3 = pick(rng);
                        if (r1 == i || r2 == i || r3 == i)continue;
                        std::vector<double> trial(ddim);
                        const int jr = static_cast<int>(u(rng) * ddim);
                        for (int d = 0; d < ddim; d++) {
                            if (u(rng) < 0.95 || d == jr)trial[d] = clamp(mp[r1][d] + 0.4 * (mp[r2][d] - mp[r3][d]));
                            else trial[d] = mp[i][d];
                        }
                        double tf = evalDM(dm, trial, ctx);
                        if (tf > mf[i]) { mp[i] = trial; mf[i] = tf; }
                    }
                    auto bit = std::max_element(mf.begin(), mf.end());
                    size_t bi = static_cast<size_t>(std::distance(mf.begin(), bit));
                    for (int k = 0; k < ddim; k++) cand[idx[k]] = mp[bi][k];
                }
            }
        }

        //  Phase 5: Joint local refinement
        void JointRefine(std::vector<double>& cand, unsigned seed = 42u) {
            const int tot = m_bench->getNumVariables();
            const int mn = 30;
            std::mt19937 rng(seed);
            std::normal_distribution<double> g(0.0, 0.015);
            std::uniform_real_distribution<double> u(0.0, 1.0);
            std::uniform_int_distribution<int> pick(0, mn - 1);
            std::vector<std::vector<double>> pop(mn, cand);
            std::vector<double> fit(mn);
            fit[0] = evalJoint(cand);
            for (int i = 1; i < mn; i++) {
                for (auto& v : pop[i])v = clamp(v + g(rng));
                fit[i] = evalJoint(pop[i]);
            }
            for (int gen = 0; gen < m_refine_gens; gen++) for (int i = 0; i < mn; i++) {
                int r1 = i, r2 = i, r3 = i, att = 0;
                while (r1 == i && ++att < 100)r1 = pick(rng);
                att = 0; while ((r2 == i || r2 == r1) && ++att < 100)r2 = pick(rng);
                att = 0; while ((r3 == i || r3 == r1 || r3 == r2) && ++att < 100)r3 = pick(rng);
                if (r1 == i || r2 == i || r3 == i)continue;
                std::vector<double> trial(tot);
                const int jr = static_cast<int>(u(rng) * tot);
                for (int d = 0; d < tot; d++) {
                    if (u(rng) < 0.98 || d == jr)trial[d] = clamp(pop[r1][d] + 0.3 * (pop[r2][d] - pop[r3][d]));
                    else trial[d] = pop[i][d];
                }
                double tf = evalJoint(trial);
                if (tf > fit[i]) { pop[i] = trial; fit[i] = tf; }
            }
            auto bit = std::max_element(fit.begin(), fit.end());
            cand = pop[static_cast<size_t>(std::distance(fit.begin(), bit))];
        }

        //  Collect final candidates (Pm + dm_elite + OSA + Cartesian party optima)
        std::vector<std::vector<double>> collectCandidates() const {
            std::vector<std::vector<double>> cands = m_mediating_pop;
            const int nd = m_bench->getNumDMs();
            for (int dm = 0; dm < nd; dm++)
                for (const auto& e : m_dm_elite[dm]) cands.push_back(e);
            for (const auto& e : m_osa) cands.push_back(e);

            // Inject Cartesian combinations of party optima explicitly
            const int tot = m_bench->getNumVariables();
            std::vector<int> cur(static_cast<size_t>(nd), 0);
            std::function<void(int)> dfs = [&](int dm) {
                if (dm == nd) {
                    std::vector<double> full(static_cast<size_t>(tot), 0.5);
                    for (int d = 0; d < nd; d++) {
                        if (m_party_optima[d].empty()) return;
                        const auto& idx = m_bench->getPartyIndices(d);
                        const auto& pos = m_party_optima[d][cur[d]];
                        for (size_t k = 0; k < idx.size() && k < pos.size(); k++)
                            full[idx[k]] = pos[k];
                    }
                    cands.push_back(full);
                    return;
                }
                for (int i = 0; i < static_cast<int>(m_party_optima[dm].size()); i++) {
                    cur[dm] = i; dfs(dm + 1);
                }
                };
            dfs(0);

            std::vector<double> fits(cands.size());
            for (size_t i = 0; i < cands.size(); i++) fits[i] = evalJoint(cands[i]);
            std::vector<size_t> ord(cands.size()); std::iota(ord.begin(), ord.end(), 0);
            std::sort(ord.begin(), ord.end(), [&](size_t a, size_t b) {return fits[a] > fits[b]; });
            std::vector<std::vector<double>> sel;
            for (size_t idx : ord) {
                bool close = false;
                for (const auto& ex : sel) {
                    double d2 = 0;
                    for (size_t d = 0; d < cands[idx].size(); d++) { double dv = cands[idx][d] - ex[d]; d2 += dv * dv; }
                    if (std::sqrt(d2) < 0.03) { close = true; break; }
                }
                if (!close) sel.push_back(cands[idx]);
            }
            return sel;
        }

        // Visualization helpers
        void savePartyLandscape(int dm, const std::vector<std::vector<double>>& ctx,
            const std::string& path) const {
            const int res = 100;
            const auto& idx = m_bench->getPartyIndices(dm);
            if (idx.size() < 2) return;
            std::ofstream f(path);
            f << std::fixed << std::setprecision(6);
            f << "# Party " << dm << " landscape (vars " << idx[0] << "," << idx[1] << ")\n";
            for (int i = 0; i < res; i++) {
                for (int j = 0; j < res; j++) {
                    std::vector<double> my_vars(m_bench->getPartyDim(dm), 0.5);
                    my_vars[0] = static_cast<double>(i) / (res - 1);
                    my_vars[1] = static_cast<double>(j) / (res - 1);
                    double fit = evalDM(dm, my_vars, ctx);
                    f << my_vars[0] << " " << my_vars[1] << " " << fit << "\n";
                }
            }
        }

        void saveFullLandscape(const std::string& out_dir) const {
            std::filesystem::create_directories(out_dir);
            const int res = 100;
            const int dim = m_bench->getNumVariables();
            for (int d1 = 0; d1 < dim; d1++) for (int d2 = d1 + 1; d2 < dim; d2++) {
                std::ostringstream fn;
                fn << out_dir << "/slice_d" << d1 << "_d" << d2 << ".txt";
                std::ofstream f(fn.str());
                f << std::fixed << std::setprecision(6);
                f << "# x" << d1 << " x" << d2 << " f [others=0.5]\n";
                for (int i = 0; i < res; i++) for (int j = 0; j < res; j++) {
                    std::vector<double> v(dim, 0.5);
                    v[d1] = static_cast<double>(i) / (res - 1);
                    v[d2] = static_cast<double>(j) / (res - 1);
                    double fv = evalJoint(v);
                    f << v[d1] << " " << v[d2] << " " << fv << "\n";
                }
            }
        }

    public:
        // Constructor
        MPM_CoEA(std::shared_ptr<MPMMO_Benchmark> bench,
            std::shared_ptr<Environment> env,
            int    K = 10,
            int    niche_pop = 30,
            int    med_pop = 400,
            int    party_gens = 150,
            int    coevo_gens = 600,
            int    refine_gens = 120,
            double sigma = 0.15,
            double de_F = 0.5,
            double de_CR = 0.9,
            double sharing_r = 0.30,
            double sharing_alpha = 1.0,
            const std::string& out_dir = "Visualization/solutions")
            : m_bench(bench), m_env(env),
            m_K(K), m_niche_pop(niche_pop), m_med_pop_size(med_pop),
            m_party_gens(party_gens), m_coevo_gens(coevo_gens),
            m_refine_gens(refine_gens), m_sigma(sigma),
            m_de_F(de_F), m_de_CR(de_CR),
            m_sharing_radius(sharing_r), m_sharing_alpha(sharing_alpha),
            m_output_dir(out_dir),
            m_osa_radius(sigma)
        {
            std::filesystem::create_directories(out_dir);
            const int nd = m_bench->getNumDMs();
            m_party_optima.resize(nd);  m_party_opt_fit.resize(nd);
            m_competing_pops.resize(nd); m_competing_fit.resize(nd);
            m_dm_elite.resize(nd);      m_dm_elite_fit.resize(nd);
            m_context_archive.resize(nd); m_context_fit.resize(nd);
        }

        //  Main run
        void run() {
            const int nd = m_bench->getNumDMs();

            // Pre-run landscapes
            auto neutral = neutralCtx();
            for (int dm = 0; dm < nd; dm++) {
                savePartyLandscape(dm, neutral, m_output_dir + "/party" + std::to_string(dm) + "_landscape_before.txt");
            }
            saveFullLandscape(m_output_dir + "/landscape_full_before");

            // Phase 1a
            std::cout << "\n>>> [MPM-CoEA] Phase 1a: NBC initialization per party\n";
            for (int dm = 0; dm < nd; dm++) NBC_init(dm, neutral);

            // Phase 1b
            std::cout << ">>> [MPM-CoEA] Phase 1b: Per-party clearing-DE ("
                << m_party_gens << " gens)\n";
            for (int dm = 0; dm < nd; dm++) ClearingDE(dm, neutral);

            // Context crossing (keep only K best diverse sub-pops)
            for (int dm = 0; dm < nd; dm++) {
                int other = 1 - dm;
                if (other < 0 || other >= nd) continue;
                int n_ctx = static_cast<int>(m_party_optima[other].size());
                if (n_ctx == 0) continue;

                auto orig_pops = m_competing_pops[dm];
                auto orig_fits = m_competing_fit[dm];

                std::vector<std::vector<std::vector<double>>> all_pops;
                std::vector<std::vector<double>> all_fits;
                for (size_t sp = 0; sp < orig_pops.size(); ++sp) {
                    all_pops.push_back(orig_pops[sp]);
                    all_fits.push_back(orig_fits[sp]);
                }

                for (int ci = 0; ci < std::min(n_ctx, m_K); ci++) {
                    auto ctx = neutral;
                    ctx[other] = m_party_optima[other][ci];
                    NBC_init(dm, ctx);
                    ClearingDE(dm, ctx, static_cast<unsigned>(dm * 31337u + ci * 9001u));
                    for (size_t sp = 0; sp < m_competing_pops[dm].size(); ++sp) {
                        all_pops.push_back(m_competing_pops[dm][sp]);
                        all_fits.push_back(m_competing_fit[dm][sp]);
                    }
                }

                // Select best K diverse sub-pops
                std::vector<std::pair<double, size_t>> best_scores;
                for (size_t sp = 0; sp < all_pops.size(); ++sp) {
                    auto it = std::max_element(all_fits[sp].begin(), all_fits[sp].end());
                    if (*it > DEAD + 1.0) best_scores.push_back({ *it, sp });
                }
                std::sort(best_scores.begin(), best_scores.end(),
                    [](const auto& a, const auto& b) { return a.first > b.first; });

                m_competing_pops[dm].clear();
                m_competing_fit[dm].clear();
                const double sep_thresh = m_sigma * 0.5;

                for (auto& [score, sp_idx] : best_scores) {
                    if (static_cast<int>(m_competing_pops[dm].size()) >= m_K) break;
                    auto& cand_pop = all_pops[sp_idx];
                    auto& cand_fit = all_fits[sp_idx];
                    auto bit = std::max_element(cand_fit.begin(), cand_fit.end());
                    size_t bi = static_cast<size_t>(std::distance(cand_fit.begin(), bit));
                    const auto& best_sol = cand_pop[bi];

                    bool too_close = false;
                    for (size_t ex = 0; ex < m_competing_pops[dm].size(); ++ex) {
                        auto& ex_fit = m_competing_fit[dm][ex];
                        auto ex_it = std::max_element(ex_fit.begin(), ex_fit.end());
                        size_t ex_bi = static_cast<size_t>(std::distance(ex_fit.begin(), ex_it));
                        double d2 = 0.0;
                        for (int d = 0; d < m_bench->getPartyDim(dm); ++d) {
                            double dv = best_sol[d] - m_competing_pops[dm][ex][ex_bi][d];
                            d2 += dv * dv;
                        }
                        if (std::sqrt(d2) < sep_thresh) { too_close = true; break; }
                    }
                    if (!too_close) {
                        m_competing_pops[dm].push_back(cand_pop);
                        m_competing_fit[dm].push_back(cand_fit);
                    }
                }

                m_party_optima[dm].clear();
                m_party_opt_fit[dm].clear();
                for (size_t sp = 0; sp < m_competing_pops[dm].size(); ++sp) {
                    auto& fit = m_competing_fit[dm][sp];
                    auto it = std::max_element(fit.begin(), fit.end());
                    if (*it > DEAD + 1.0) {
                        size_t bi = static_cast<size_t>(std::distance(fit.begin(), it));
                        m_party_optima[dm].push_back(m_competing_pops[dm][sp][bi]);
                        m_party_opt_fit[dm].push_back(*it);
                    }
                }

                std::cout << "  Party " << dm << ": " << m_party_optima[dm].size()
                    << " optima after context crossing (selected from " << all_pops.size() << ")\n";
            }

            std::cout << "  Per-party optima:\n";
            for (int dm = 0; dm < nd; dm++)
                for (size_t k = 0; k < m_party_optima[dm].size(); k++)
                    std::cout << "    [P" << dm << " opt" << k << "] f="
                    << std::fixed << std::setprecision(4) << m_party_opt_fit[dm][k] << "\n";

            // Phase 2
            std::cout << "\n>>> [MPM-CoEA] Phase 2: Cartesian assembly :: seed mediating pop\n";
            CartesianSeed();

            // Phase 3: Co-evolutionary loop
            std::cout << "\n>>> [MPM-CoEA] Phase 3: Co-evolutionary loop ("
                << m_coevo_gens << " gens, logged every 5)\n";
            m_best_prev = *std::max_element(m_mediating_fit.begin(), m_mediating_fit.end());
            m_stagnation_cnt = 0;

            m_conv_log.push_back({ 0, m_best_prev });
            std::cout << "  Gen    0 | best=" << std::fixed << std::setprecision(4) << m_best_prev
                << " | stag=0 | restarts=0"
                << " | osa_size=" << m_osa.size() << "\n";

            for (m_generation = 1; m_generation <= m_coevo_gens; m_generation++) {
                adaptOsaRadius();
                evolveCompeting();
                updateContextArchives();
                evolveMediating();
                applyDMPreservation();
                checkAndRestart();

                if (m_generation % 10 == 0) {
                    saveSolutionsAtGeneration(m_generation);
                }

                if (m_generation % 5 == 0) {
                    double best = *std::max_element(m_mediating_fit.begin(), m_mediating_fit.end());
                    std::cout << "  Gen " << std::setw(4) << m_generation
                        << " | best=" << std::fixed << std::setprecision(4) << best
                        << " | stag=" << m_stagnation_cnt
                        << " | restarts=" << m_restart_count
                        << " | osa_size=" << m_osa.size() << "\n";
                }
            }

            // Phase 4 + 5: Refinement
            std::cout << "\n>>> [MPM-CoEA] Phase 4-5: Refinement\n";
            auto finals = collectCandidates();
            finals.resize(std::min(static_cast<int>(finals.size()), 120));
            for (size_t i = 0; i < finals.size(); i++) ContextRefine(finals[i], 3);
            for (size_t i = 0; i < finals.size(); i++) {
                JointRefine(finals[i], static_cast<unsigned>(i * 13u + 42u));
                osaInsert(finals[i], evalJoint(finals[i]));
            }
            for (size_t i = 0; i < finals.size() && i < m_mediating_pop.size(); i++) {
                m_mediating_pop[i] = finals[i];
                m_mediating_fit[i] = evalJoint(finals[i]);
            }
            double best = *std::max_element(m_mediating_fit.begin(), m_mediating_fit.end());
            std::cout << "  Best after refinement: " << std::fixed << std::setprecision(4) << best << "\n";
            std::cout << "  OSA archive size: " << m_osa.size() << "\n";

            // Post-run landscapes
            for (int dm = 0; dm < nd; dm++) {
                auto best_ctx = ctxFromJoint(dm, getBestSolution());
                savePartyLandscape(dm, best_ctx, m_output_dir + "/party" + std::to_string(dm) + "_landscape_after.txt");
            }
            saveFullLandscape(m_output_dir + "/landscape_full_after");
        }

        // Accessors
        double getBestFitness() const {
            if (m_mediating_fit.empty()) return DEAD;
            return *std::max_element(m_mediating_fit.begin(), m_mediating_fit.end());
        }

        std::vector<double> getBestSolution() const {
            if (m_mediating_pop.empty()) return {};
            auto it = std::max_element(m_mediating_fit.begin(), m_mediating_fit.end());
            return m_mediating_pop[static_cast<size_t>(std::distance(m_mediating_fit.begin(), it))];
        }

        std::vector<std::vector<double>> getAllCandidates() const {
            return collectCandidates();
        }

        void saveConvergenceCurve() const {
            std::ofstream out(m_output_dir + "/convergence.txt");
            if (!out.is_open()) return;
            out << "# generation best_joint_fitness\n";
            for (const auto& p : m_conv_log)
                out << p.first << " " << std::fixed << std::setprecision(6) << p.second << "\n";
        }

        void saveSolutionsAtGeneration(int tag) const {
            std::string fname = m_output_dir + "/gen_" + std::to_string(tag) + "_solutions.txt";
            std::ofstream out(fname);
            if (!out.is_open()) return;
            out << std::fixed << std::setprecision(6) << "# gen=" << tag << " x0..xd f\n";
            for (size_t i = 0; i < m_mediating_pop.size(); i++) {
                for (double v : m_mediating_pop[i]) out << v << " ";
                out << m_mediating_fit[i] << "\n";
            }
            out << "# OSA solutions\n";
            for (size_t i = 0; i < m_osa.size(); i++) {
                for (double v : m_osa[i]) out << v << " ";
                out << m_osa_fit[i] << "\n";
            }
        }
    };

    //  Benchmark creation
    void create4D2PartyProblem(const std::string& problem_name = "complex_4d_2p") {
        using namespace free_peaks;
        const std::string pd = "multiparty_multimodal/", sd = "sop/";
        for (auto& p : { pd, sd,
                        std::string("subproblem/") + sd,
                        std::string("subproblem/function/one_peak/") + sd })
            std::filesystem::create_directories(
                g_working_directory + "instance/problem/continuous/free_peaks/" + p);

        const int numDim = 4, numObj = 1, numCon = 0;
        std::shared_ptr<ofec::Random> rnd(new ofec::Random(0.123));
        FreePeaks::registerFP();
        std::shared_ptr<Environment> env(generateEnvironmentByFactory("free_peaks"));
        env->recordInputParameters(); env->initialize();
        env->setProblem(generateProblemByFactory("free_peaks"));
        auto* fp = dynamic_cast<FreePeaks*>(env->problem());

        ParameterMap pm;
        pm["generation_type"] = std::string("assigned");
        pm["dataFile1"] = pd + problem_name + ".txt";
        fp->inputParameters().input(pm); fp->initialize(env.get());
        fp->setRandom(rnd); fp->setSizes(numDim, numObj, numCon);

        fp->setKDtree({ {"root",{
            {"peak1",0.20},{"peak2",0.16},{"peak3",0.13},{"peak4",0.12},
            {"peak5",0.11},{"peak6",0.10},{"peak7",0.10},{"peak8",0.08}}} });

        struct PeakSpec {
            std::vector<double> center;
            double h;
            std::string shape;
            double cond;
            int party;
            double mag, rot;
        };

        const std::vector < PeakSpec > peaks = {
            {{  0,  0,  0,  0}, 90.0,"s1",   1.0, 0, 0.0, 0.0},
            {{ 40,-30, 35,-25}, 75.0,"s2",  10.0, 0, 2.0, 0.2},
            {{-50, 40,-45, 35}, 65.0,"s7",  10.0, 1, 2.0, 0.2},
            {{ 25, 25, 25, 25}, 55.0,"s1",  50.0, 0, 2.0, 0.0},
            {{-30,-30,-30,-30}, 50.0,"s2",   1.0, 1, 2.0, 0.0},
            {{ 60,-50, 55,-45}, 42.0,"s2", 100.0, 1, 2.0, 0.0},
            {{-20, 60,-15, 55}, 35.0,"s1",  50.0, 0, 2.0, 0.0},
            {{ 80, 80, 75, 80}, 80.0,"s7", 100.0, 1, 2.0, 0.2},
        };
        const std::vector<std::string> sp_names = {
            "peak1","peak2","peak3","peak4","peak5","peak6","peak7","peak8" };

        for (size_t s = 0; s < sp_names.size(); s++) {
            const auto& sname = sp_names[s]; const auto& pk = peaks[s];
            ParameterMap spm;
            spm["subspace"] = sname; spm["generation_type"] = std::string("assigned");
            spm["dataFile1"] = sd + problem_name + "_" + sname + ".txt";
            auto subpro(Subproblem::create()); subpro->initialize(spm, fp);

            {
                auto dis(FactoryFP<DistanceBase>::produce("Euclidean")); ParameterMap dp;
                dis->initialize(fp, sname, dp); subpro->setDistance(dis);
            }

            {
                ParameterMap fpm; fpm["generation_type"] = std::string("assigned");
                fpm["dataFile1"] = sd + problem_name + "_" + sname + "_fn.txt";
                auto func(FactoryFP<FunctionBase>::produce("one_peak"));
                func->initialize(fp, sname, fpm);
                auto* opf = dynamic_cast<OnePeakFunction*>(func);
                auto op(FactoryFP<OnePeakBase>::produce(pk.shape));
                ParameterMap opp; opp["center_type"] = std::string("assigned");
                opp["height"] = static_cast<Real>(pk.h);
                std::vector<Real> c(pk.center.size());
                for (size_t d = 0; d < pk.center.size(); d++)
                    c[d] = static_cast<Real>(std::max(-100.0, std::min(100.0, pk.center[d])));
                opp["center_postion"] = c; op->initialize(fp, sname, opp);
                opf->addOnePeaks(op); subpro->setFunction(func);
            }

            // Transform 1: Party-bias shift + optional rotation
            {
                auto t(FactoryFP<X_TransformBase>::produce("MapXPartyBias")); ParameterMap tp;
                tp["party_id"] = pk.party; tp["magnitude"] = static_cast<Real>(pk.mag);
                tp["rotation_angle"] = static_cast<Real>(pk.rot); tp["condition_number"] = Real(1.0);
                t->initialize(fp, sname, tp); subpro->addVariableTransform(t);
            }

            // Transform 2: Ill-conditioning
            {
                auto t(FactoryFP<X_TransformBase>::produce("MapXIllConditioning")); ParameterMap tp;
                tp["condition"] = static_cast<Real>(pk.cond); t->initialize(fp, sname, tp);
                subpro->addVariableTransform(t);
            }

            // Transform 3: Irregularity
            {
                auto t(FactoryFP<X_TransformBase>::produce("MapXIrregularity")); ParameterMap tp;
                t->initialize(fp, sname, tp); subpro->addVariableTransform(t);
            }

            // Transform 4: Linkage
            {
                auto t(FactoryFP<X_TransformBase>::produce("MapXLinkage")); ParameterMap tp;
                tp["beta"] = Real(0.05); t->initialize(fp, sname, tp); subpro->addVariableTransform(t);
            }

            // Transform 5: Full random rotation
            {
                auto t(FactoryFP<X_TransformBase>::produce("MapXRotation")); ParameterMap tp;
                t->initialize(fp, sname, tp); subpro->addVariableTransform(t);
            }

            // Transform 6: Asymmetry
            {
                auto t(FactoryFP<X_TransformBase>::produce("MapXAsymmetry"));
                ParameterMap tp;
                tp["beta"] = Real(0.5);
                t->initialize(fp, sname, tp);
                subpro->addVariableTransform(t);
            }

            {
                auto ot(FactoryFP<Y_TransformBase>::produce("map_objective")); ParameterMap otp;
                ot->initialize(fp, sname, otp); subpro->addObjectiveTransform(ot);
            }
            fp->setSubproblem(sname, subpro);
        }
        fp->bindData();
        fp->inputParameters().at("generation_type")->setValue("read_file");
        for (auto& it : fp->subspaceTree().name_box_subproblem) if (it.second.second) {
            it.second.second->inputParameters().at("generation_type")->setValue("read_file");
            it.second.second->function()->inputParameters().at("generation_type")->setValue("read_file");
        }
        fp->recordInputParameters(); fp->outputTotalFile();
        std::cout << ">>> Problem '" << problem_name << "' created: 8 peaks, dim=4, 2 parties\n";
    }

    //  Helper functions
    void inspectPeakCenters(FreePeaks* problem, Environment* env) {
        const auto* opt = problem->optima();
        if (!opt || opt->numberSolutions() == 0) { std::cout << "[INSPECT] No optima.\n"; return; }
        std::cout << "\n>>> MANUAL PEAK CENTER INSPECTION\n    Total peaks: " << opt->numberSolutions() << "\n";
        double mx = std::numeric_limits<double>::lowest(); int bi = -1;
        for (int i = 0; i < static_cast<int>(opt->numberSolutions()); i++) {
            const auto& sol = opt->solution(i); const auto& v = sol.variable().vector();
            std::cout << "    Peak " << (i + 1) << ": [";
            for (size_t d = 0; d < v.size(); d++)
                std::cout << std::fixed << std::setprecision(4)
                << static_cast<double>(v[d])
                << (d + 1 < v.size() ? ", " : "");
            std::cout << "] -> " << std::fixed << std::setprecision(4) << sol.objective(0) << "\n";
            if (sol.objective(0) > mx) { mx = sol.objective(0); bi = i; }
        }
        std::cout << "    Global optimum: Peak " << (bi + 1) << " with fitness "
            << std::fixed << std::setprecision(4) << mx << "\n";
        std::cout << "\n    Direct evaluation verification:\n";
        for (int i = 0; i < static_cast<int>(opt->numberSolutions()); i++) {
            const auto& sol = opt->solution(i); const auto& v = sol.variable().vector();
            std::unique_ptr<SolutionBase> ts(problem->createSolution());
            auto& sv = dynamic_cast<Solution<VariableVector<Real>>&>(*ts).variable().vector();
            for (size_t d = 0; d < v.size() && d < sv.size(); d++) sv[d] = v[d];
            ts->evaluate(env, false);
            double ef = ts->objective(0);
            std::cout << "    Peak " << (i + 1) << " direct eval: " << std::fixed << std::setprecision(4) << ef;
            if (std::abs(ef - sol.objective(0)) > 0.1) std::cout << " [WARNING stored=" << sol.objective(0) << "]";
            std::cout << "\n";
        }
    }

    void outputLandscapeSlices(Environment* env, const std::string& out_dir) {
        std::filesystem::create_directories(out_dir);
        const int res = 200;
        auto* cp = CAST_CONOP(env->problem()); const int dim = static_cast<int>(cp->numberVariables());
        std::cout << ">>> Sampling landscape slices (res=" << res << ")..." << std::flush;
        for (int d1 = 0; d1 < dim; d1++) for (int d2 = d1 + 1; d2 < dim; d2++) {
            std::ostringstream fn; fn << out_dir << "/slice_d" << d1 << "_d" << d2 << ".txt";
            std::ofstream f(fn.str()); f << std::fixed << std::setprecision(6);
            f << "# x" << d1 << " x" << d2 << " f [others=0.5]\n";
            for (int i = 0; i < res; i++) for (int j = 0; j < res; j++) {
                auto sol = env->problem()->createSolution();
                auto& v = dynamic_cast<Solution<VariableVector<Real>>&>(*sol).variable().vector();
                v.assign(dim, Real(0.5)); v[d1] = Real(i) / Real(res - 1); v[d2] = Real(j) / Real(res - 1);
                sol->evaluate(env, false);
                f << (double(i) / (res - 1)) << " " << (double(j) / (res - 1)) << " " << sol->objective(0) << "\n";
            }
        }
        std::cout << " done.\n";
    }

    //  Main experiment
    void runMPMCoEAExperiment() {
        ofec::registerInstance();
        ofec::g_working_directory =
            R"(E:\HITSZ\Research\Multimodal_Multiparty_Optimization\ThesisProject\Data\ofec_data_new\)";
        for (auto& c : ofec::g_working_directory) if (c == '\\')c = '/';
        const std::string vis_base =
            R"(E:\HITSZ\Research\Multimodal_Multiparty_Optimization\ThesisProject\Visualization)";

        const std::string problem_name = "complex_4d_2p";
        std::cout << "\n>>> [MPM-CoEA] Creating 4D 2-party benchmark...\n";
        create4D2PartyProblem(problem_name);

        auto* env_raw = generateEnvironmentByFactory("free_peaks");
        env_raw->recordInputParameters(); env_raw->initialize();
        env_raw->setProblem(generateProblemByFactory("free_peaks"));
        ParameterMap lpm;
        lpm["generation_type"] = std::string("read_file");
        lpm["dataFile1"] = std::string("multiparty_multimodal/") + problem_name + ".txt";
        env_raw->problem()->inputParameters().input(lpm);
        env_raw->problem()->recordInputParameters();
        env_raw->initializeProblem(0.5);

        auto* fp = CAST_FPs(env_raw->problem());
        if (!fp) throw std::runtime_error("Problem is not FreePeaks.");
        const int nopt = fp->optima() ? static_cast<int>(fp->optima()->numberSolutions()) : 0;
        std::cout << ">>> Loaded: dim=" << fp->numberVariables() << " | parties=2 | num_optima=" << nopt << "\n";
        inspectPeakCenters(fp, env_raw);
        outputLandscapeSlices(env_raw, vis_base + "/landscape_before");

        auto bench = std::make_shared<MPMMO_Benchmark>(fp, 2, env_raw);
        auto env_sp = std::shared_ptr<Environment>(env_raw, [](Environment*) {});

        MPM_CoEA solver(bench, env_sp,
            /*K*/           10,
            /*niche_pop*/   30,
            /*med_pop*/     400,
            /*party_gens*/  150,
            /*coevo_gens*/  600,
            /*refine_gens*/ 120,
            /*sigma*/       0.15,
            /*de_F*/        0.5,
            /*de_CR*/       0.9,
            /*sharing_r*/   0.30,
            /*sharing_a*/   1.0,
            vis_base + "/solutions");

        solver.run();
        solver.saveSolutionsAtGeneration(600);
        solver.saveConvergenceCurve();
        outputLandscapeSlices(env_raw, vis_base + "/landscape_after");

        auto candidates = solver.getAllCandidates();

        std::cout << "\n>>> CEC-2015 Metrics (STRICT: eps_x=0.01, eps_f=1.0)\n";
        auto ms = evaluateCEC2015Metrics(candidates, fp, env_raw, 0.01, 1.0);
        std::cout << "    SR  : " << ms.SR << "%\n    ANOF: " << ms.ANOF
            << " (" << static_cast<int>(std::round(ms.ANOF * nopt)) << "/" << nopt << " found)\n"
            << "    MPR : " << ms.MPR << "\n    SP  : " << ms.SP << "\n";

        std::cout << "\n>>> CEC-2015 Metrics (LOOSE: eps_x=0.05, eps_f=2.0)\n";
        auto ml = evaluateCEC2015Metrics(candidates, fp, env_raw, 0.05, 2.0);
        std::cout << "    SR  : " << ml.SR << "%\n    ANOF: " << ml.ANOF
            << " (" << static_cast<int>(std::round(ml.ANOF * nopt)) << "/" << nopt << " found)\n"
            << "    MPR : " << ml.MPR << "\n    SP  : " << ml.SP << "\n";

        auto best = solver.getBestSolution();
        std::cout << ">>> Best consensus solution: [";
        for (size_t i = 0; i < best.size(); i++)
            std::cout << std::fixed << std::setprecision(4) << best[i] << (i + 1 < best.size() ? ", " : "");
        std::cout << "]  f = " << std::fixed << std::setprecision(4) << solver.getBestFitness() << "\n";
        std::cout << "\n>>> [MPM-CoEA] Complete.\n";
    }

} // namespace ofec
#endif // OFEC_MPMCOEA_SOLVER_HPP