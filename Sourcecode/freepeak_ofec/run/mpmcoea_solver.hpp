#ifndef OFEC_MPMCOEA_SOLVER_HPP
#define OFEC_MPMCOEA_SOLVER_HPP

#include "../core/global.h"
#include "interface.h"
#include "../instance/problem/continuous/free_peaks/free_peaks.h"
#include "../instance/problem/continuous/free_peaks/free_peaks_multiparty.h"
#include "../instance/problem/continuous/free_peaks/subproblem/subproblem.h"
#include "../instance/problem/continuous/free_peaks/subproblem/function/one_peak_function.h"
#include "../core/environment/environment.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace ofec {

    //  Nearest-Better Clustering (NBC)
    inline std::vector<std::vector<size_t>> nearestBetterClustering(
        const std::vector<std::vector<double>>& pop,
        const std::vector<double>& fit,
        double phi)
    {
        const int n = static_cast<int>(pop.size());
        if (n == 0) return {};

        std::vector<int>    nbi(n, -1);
        std::vector<double> nbi_dist(n, std::numeric_limits<double>::max());
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (fit[j] <= fit[i] + 1e-12) continue;
                double d2 = 0.0;
                for (size_t k = 0; k < pop[i].size(); ++k) {
                    double dv = pop[i][k] - pop[j][k]; d2 += dv * dv;
                }
                double d = std::sqrt(d2);
                if (d < nbi_dist[i]) { nbi_dist[i] = d; nbi[i] = j; }
            }
        }

        double mean_len = 0.0; int cnt = 0;
        for (int i = 0; i < n; ++i)
            if (nbi[i] >= 0) { mean_len += nbi_dist[i]; ++cnt; }
        if (cnt > 0) mean_len /= cnt;
        const double thresh = phi * mean_len;

        std::vector<std::vector<int>> adj(n);
        for (int i = 0; i < n; ++i) {
            if (nbi[i] >= 0 && nbi_dist[i] <= thresh) {
                adj[i].push_back(nbi[i]);
                adj[nbi[i]].push_back(i);
            }
        }

        std::vector<int> visited(n, 0);
        std::vector<std::vector<size_t>> clusters;
        for (int s = 0; s < n; ++s) {
            if (visited[s]) continue;
            std::vector<size_t> cluster;
            std::vector<int> q = { s };
            visited[s] = 1;
            while (!q.empty()) {
                int u = q.back(); q.pop_back();
                cluster.push_back(static_cast<size_t>(u));
                for (int v : adj[u]) if (!visited[v]) { visited[v] = 1; q.push_back(v); }
            }
            clusters.push_back(std::move(cluster));
        }
        return clusters;
    }

    //  MPMMO_Benchmark - wraps FreePeaksMultiParty
    class MPMMO_Benchmark {
    public:
        struct Eval {
            std::array<double, 2> obj{ 0.0, 0.0 };
            double consensus = 0.0;
        };

        MPMMO_Benchmark(FreePeaksMultiParty* problem, Environment* env)
            : m_problem(problem), m_env(env)
            , m_dim(static_cast<int>(problem->numberVariables()))
        {
        }

        Eval evaluate(const std::vector<double>& x) const {
            auto sol = std::unique_ptr<SolutionBase>(m_problem->createSolution());
            auto& v = dynamic_cast<Solution<VariableVector<Real>>&>(*sol)
                .variable().vector();
            for (int d = 0; d < m_dim; ++d)
                v[d] = Real(std::max(0.0, std::min(1.0, x[d])));
            sol->evaluate(m_env, false);
            Eval r;
            r.obj[0] = static_cast<double>(sol->objective(0));
            r.obj[1] = static_cast<double>(sol->objective(1));
            r.consensus = std::min(r.obj[0], r.obj[1]);
            return r;
        }

        double evaluateParty(int p, const std::vector<double>& x) const {
            return evaluate(x).obj[static_cast<size_t>(p)];
        }

        //  Dimension-normalised party evaluation for Phase 1b search.
        //
        //  Gaussian peaks exp(-r^2 / 2*sig^2) become exponentially narrow in
        //  high D: at a fixed distance from the peak centre, fitness drops as
        //  exp(-D*delta^2 / 2*sig^2). For D=10 the peaks are invisible to DE.
        //
        //  Rescaling:  f_adj = 100 * (f/100)^(2/D)
        //    - At D=2: exponent=1.0, no change.
        //    - At D=10: exponent=0.2, stretches the gradient near-zero values
        //      up so Phase 1b clearing-DE can detect peak structure.
        //  Used ONLY for per-party search (NBC init + ClearingDE).
        //  Archives and metrics use raw evaluate() values.
        double evaluatePartyScaled(int p, const std::vector<double>& x) const {
            double f = evaluateParty(p, x);
            if (m_dim <= 2 || f < 1e-12) return f;
            const double scale = 2.0 / static_cast<double>(m_dim);
            return 100.0 * std::pow(f / 100.0, scale);
        }

        double consensusFitness(const std::vector<double>& x) const {
            return evaluate(x).consensus;
        }

        int numParties()      const { return 2; }
        int dim()             const { return m_dim; }

        FreePeaksMultiParty* problem() const { return m_problem; }
        Environment* env()     const { return m_env; }

        const std::vector<std::vector<Real>>& sharedOptima() const {
            return m_problem->currentSpec().shared_optima;
        }
        std::string problemName() const { return m_problem->currentSpec().name; }
        std::string feature()     const { return m_problem->currentSpec().feature; }

    private:
        FreePeaksMultiParty* m_problem;
        Environment* m_env;
        int                  m_dim;
    };

    //  CEC-2015 Metrics
    struct CEC2015Result { double SR, ANOF, MPR; };

    inline CEC2015Result evaluateCEC2015Metrics(
        const std::vector<std::vector<double>>& candidates,
        const MPMMO_Benchmark& bench,
        double eps_x,
        double eps_f)
    {
        const auto& shared_opts = bench.sharedOptima();
        const int nopt = static_cast<int>(shared_opts.size());
        if (nopt == 0) return { 0.0, 0.0, 0.0 };

        std::vector<double> opt_c(nopt);
        for (int k = 0; k < nopt; ++k) {
            std::vector<double> xd(shared_opts[k].begin(), shared_opts[k].end());
            opt_c[k] = bench.consensusFitness(xd);
        }

        std::vector<bool> covered(nopt, false);
        for (const auto& c : candidates) {
            double c_c = bench.consensusFitness(c);
            for (int k = 0; k < nopt; ++k) {
                if (covered[k]) continue;
                double d2 = 0.0;
                for (size_t d = 0; d < c.size() && d < shared_opts[k].size(); ++d) {
                    double dv = c[d] - static_cast<double>(shared_opts[k][d]);
                    d2 += dv * dv;
                }
                if (std::sqrt(d2) < eps_x && c_c >= opt_c[k] - eps_f)
                    covered[k] = true;
            }
        }

        int n_found = static_cast<int>(std::count(covered.begin(), covered.end(), true));
        return {
            (n_found == nopt) ? 100.0 : 0.0,
            static_cast<double>(n_found),
            static_cast<double>(n_found) / static_cast<double>(nopt)
        };
    }

    //  MPM-CoEA
    class MPM_CoEA {
    private:
        std::shared_ptr<MPMMO_Benchmark> m_bench;
        std::shared_ptr<Environment>     m_env_keeper;

        int    m_K;
        int    m_niche_pop;
        int    m_med_pop_size;
        int    m_party_gens;
        int    m_coevo_gens;
        int    m_refine_gens;
        double m_sigma;
        double m_de_F;
        double m_de_CR;
        double m_sharing_alpha;

        static constexpr int    STAGNATION_LIMIT = 60;
        static constexpr double DEAD = -1e30;
        static constexpr double CLEARING_RADIUS = 0.12;

        // OSA
        static constexpr double OSA_MIN_QUALITY = 25.0;
        static constexpr int    OSA_MAX = 120;

        // Consensus Archive constants
        static constexpr double CA_RADIUS = 0.03;
        static constexpr int    CA_MAX = 200;

        // dim-dependent CA quality threshold (set in constructor)
        double m_ca_min_quality;

        // Landscape vis resolution
        static constexpr int    VIS_RES = 100;

        // Per-party competing subpopulations
        std::vector<std::vector<std::vector<std::vector<double>>>> m_competing_pops;
        std::vector<std::vector<std::vector<double>>>              m_competing_fit;

        // Best solutions per party
        std::vector<std::vector<std::vector<double>>> m_party_optima;
        std::vector<std::vector<double>>              m_party_opt_fit;

        // DM elites
        std::vector<std::vector<std::vector<double>>> m_dm_elite;
        std::vector<std::vector<double>>              m_dm_elite_fit;

        // Mediating population
        std::vector<std::vector<double>> m_mediating_pop;
        std::vector<double>              m_mediating_fit;
        std::vector<std::vector<double>> m_mediating_obj;

        // OSA
        std::vector<std::vector<double>> m_osa;
        std::vector<double>              m_osa_fit;
        double                           m_osa_radius;

        // Consensus Archive
        std::vector<std::vector<double>> m_ca;
        std::vector<double>              m_ca_fit;

        // Known basins from Phase 1b party optima
        // m_known_basins[b] = basin centre (best optima from any party near that peak)
        // m_basin_ca_slot[b] = index into m_ca of this basin's representative (-1 = uncovered)
        std::vector<std::vector<double>> m_known_basins;
        std::vector<int>                 m_basin_ca_slot;

        std::vector<std::pair<int, double>> m_conv_log;

        int    m_generation = 0;
        int    m_restart_count = 0;
        int    m_stagnation_cnt = 0;
        double m_best_prev = DEAD;
        std::string m_out_dir;

        //  Helpers
        static double clamp01(double v) { return std::max(0.0, std::min(1.0, v)); }

        MPMMO_Benchmark::Eval doEval(const std::vector<double>& x) const {
            return m_bench->evaluate(x);
        }
        double evalParty(const std::vector<double>& x, int p) const {
            return doEval(x).obj[static_cast<size_t>(p)];
        }
        // scaled party eval for Phase 1b search
        double evalPartyScaled(const std::vector<double>& x, int p) const {
            return m_bench->evaluatePartyScaled(p, x);
        }
        double evalConsensus(const std::vector<double>& x) const {
            return doEval(x).consensus;
        }

        // Return Euclidean distance between two vectors
        static double euclidDist(const std::vector<double>& a, const std::vector<double>& b) {
            double d2 = 0.0;
            for (size_t i = 0; i < a.size(); ++i) {
                double dv = a[i] - b[i]; d2 += dv * dv;
            }
            return std::sqrt(d2);
        }

        // index of nearest known basin to sol (-1 if no basins defined)
        int nearestBasin(const std::vector<double>& sol) const {
            if (m_known_basins.empty()) return -1;
            int best = 0;
            double best_d = euclidDist(sol, m_known_basins[0]);
            for (int b = 1; b < static_cast<int>(m_known_basins.size()); ++b) {
                double d = euclidDist(sol, m_known_basins[b]);
                if (d < best_d) { best_d = d; best = b; }
            }
            return best;
        }

        // Build known-basins list from Phase 1b party optima
        // Called once after all parties' ClearingDE is done
        void buildKnownBasins() {
            const int D = m_bench->dim();
            m_known_basins.clear();

            for (int p = 0; p < m_bench->numParties(); ++p) {
                for (const auto& opt : m_party_optima[p]) {
                    // Only add if far enough from existing basins
                    bool too_close = false;
                    for (const auto& b : m_known_basins) {
                        if (euclidDist(opt, b) < CA_RADIUS * 2.0) { too_close = true; break; }
                    }
                    if (!too_close) m_known_basins.push_back(opt);
                }
            }

            // Initialise all basin slots as uncovered
            m_basin_ca_slot.assign(m_known_basins.size(), -1);

            std::cout << "  [Basin] " << m_known_basins.size()
                << " distinct basins from Phase 1b party optima\n";
        }

        // Basin-tagged Consensus Archive insertion
        void caInsert(const std::vector<double>& sol, double c) {
            if (!std::isfinite(c) || c < m_ca_min_quality) return;

            // Find nearest CA member
            int nearest = -1;
            double min_d = std::numeric_limits<double>::max();
            for (int i = 0; i < static_cast<int>(m_ca.size()); ++i) {
                double d = euclidDist(sol, m_ca[i]);
                if (d < min_d) { min_d = d; nearest = i; }
            }

            // Close to existing: update in-place if better
            if (nearest >= 0 && min_d < CA_RADIUS) {
                if (c > m_ca_fit[nearest]) {
                    m_ca[nearest] = sol;
                    m_ca_fit[nearest] = c;
                    // Re-tag basin for this slot
                    int b = nearestBasin(sol);
                    if (b >= 0) m_basin_ca_slot[b] = nearest;
                }
                return;
            }

            // New distinct position
            int b = nearestBasin(sol);
            bool uncovered = (b >= 0 && m_basin_ca_slot[b] < 0);

            if (static_cast<int>(m_ca.size()) < CA_MAX) {
                int idx = static_cast<int>(m_ca.size());
                m_ca.push_back(sol);
                m_ca_fit.push_back(c);
                if (uncovered) m_basin_ca_slot[b] = idx;
                return;
            }

            // Archive full: find best eviction candidate.
            // Prefer to evict a redundant slot (its basin already has another slot).
            int evict = -1;
            double worst_c = std::numeric_limits<double>::max();

            for (int i = 0; i < static_cast<int>(m_ca.size()); ++i) {
                int ib = nearestBasin(m_ca[i]);
                // Is slot i the designated representative of its basin?
                bool is_sole_rep = (ib >= 0 && m_basin_ca_slot[ib] == i);
                if (is_sole_rep) continue;  // protect basin representative
                if (m_ca_fit[i] < worst_c) { worst_c = m_ca_fit[i]; evict = i; }
            }

            // If no redundant slot found (all are basin reps), fall back to worst overall
            if (evict < 0) {
                auto it = std::min_element(m_ca_fit.begin(), m_ca_fit.end());
                worst_c = *it;
                evict = static_cast<int>(std::distance(m_ca_fit.begin(), it));
            }

            // Evict if new solution is better OR if it covers an uncovered basin
            if (c > worst_c || uncovered) {
                // Fix basin slot for the evicted entry
                int old_b = nearestBasin(m_ca[evict]);
                if (old_b >= 0 && m_basin_ca_slot[old_b] == evict)
                    m_basin_ca_slot[old_b] = -1;

                m_ca[evict] = sol;
                m_ca_fit[evict] = c;
                if (uncovered) m_basin_ca_slot[b] = evict;
            }
        }

        //  OSA insertion
        void osaInsert(const std::vector<double>& sol, double c) {
            if (c <= DEAD + 1.0 || !std::isfinite(c) || c < OSA_MIN_QUALITY) return;
            double min_d = std::numeric_limits<double>::max();
            size_t close_i = 0;
            for (size_t i = 0; i < m_osa.size(); ++i) {
                double d = euclidDist(sol, m_osa[i]);
                if (d < min_d) { min_d = d; close_i = i; }
            }
            if (static_cast<int>(m_osa.size()) < OSA_MAX) {
                if (min_d >= m_osa_radius) {
                    m_osa.push_back(sol); m_osa_fit.push_back(c);
                }
                else if (!m_osa.empty() && c > m_osa_fit[close_i]) {
                    m_osa[close_i] = sol; m_osa_fit[close_i] = c;
                }
            }
            else {
                if (min_d > m_osa_radius) {
                    auto wt = std::min_element(m_osa_fit.begin(), m_osa_fit.end());
                    if (c > *wt) {
                        size_t wi = static_cast<size_t>(
                            std::distance(m_osa_fit.begin(), wt));
                        m_osa[wi] = sol; m_osa_fit[wi] = c;
                    }
                }
                else if (!m_osa.empty() && c > m_osa_fit[close_i]) {
                    m_osa[close_i] = sol; m_osa_fit[close_i] = c;
                }
            }
        }

        //  Phase 1a: NBC init per party (full D-space)
        //  Uses evalPartyScaled for fitness in NBC clustering.
        void NBC_init(int party, unsigned seed = 0) {
            const int D = m_bench->dim();
            const int N = m_K * m_niche_pop;
            std::mt19937 rng(seed ? seed
                : static_cast<unsigned>(party * 17417u + 100003u));
            std::uniform_real_distribution<double> u(0.0, 1.0);

            std::vector<std::vector<double>> pop(N, std::vector<double>(D));
            std::vector<double> fit(N);
            for (int i = 0; i < N; ++i) {
                for (auto& v : pop[i]) v = u(rng);
                // use scaled fitness for clustering so high-D peaks are detectable
                fit[i] = evalPartyScaled(pop[i], party);
            }

            auto clusters = nearestBetterClustering(pop, fit, m_sigma);

            std::sort(clusters.begin(), clusters.end(),
                [&](const std::vector<size_t>& a, const std::vector<size_t>& b) {
                    auto ba = *std::max_element(a.begin(), a.end(),
                        [&](size_t x, size_t y) { return fit[x] < fit[y]; });
                    auto bb = *std::max_element(b.begin(), b.end(),
                        [&](size_t x, size_t y) { return fit[x] < fit[y]; });
                    return fit[ba] > fit[bb];
                });
            if (static_cast<int>(clusters.size()) > m_K)
                clusters.resize(static_cast<size_t>(m_K));

            m_competing_pops[party].clear();
            m_competing_fit[party].clear();
            m_party_optima[party].clear();
            m_party_opt_fit[party].clear();

            for (const auto& cl : clusters) {
                std::vector<std::vector<double>> sp;
                std::vector<double>              sf;
                for (size_t idx : cl) { sp.push_back(pop[idx]); sf.push_back(fit[idx]); }
                auto bfit = std::max_element(sf.begin(), sf.end());
                if (*bfit > DEAD + 1.0) {
                    size_t bi = static_cast<size_t>(std::distance(sf.begin(), bfit));
                    m_party_optima[party].push_back(sp[bi]);
                    m_party_opt_fit[party].push_back(evalParty(sp[bi], party)); // raw fit for logging
                }
                m_competing_pops[party].push_back(sp);
                m_competing_fit[party].push_back(sf);
            }
        }

        //  Phase 1b: ClearingDE per party (full D-space)
        //  Uses evalPartyScaled so high-D Gaussian peaks are detectable.
        void ClearingDE(int party, unsigned seed = 0) {
            const int D = m_bench->dim();
            std::mt19937 rng(seed ? seed
                : static_cast<unsigned>(party * 31337u + 77777u));
            std::uniform_real_distribution<double> u(0.0, 1.0);
            constexpr double PARTY_CR = 0.50;

            for (int gen = 0; gen < m_party_gens; ++gen) {
                for (size_t sp = 0; sp < m_competing_pops[party].size(); ++sp) {
                    auto& pop = m_competing_pops[party][sp];
                    auto& fit = m_competing_fit[party][sp];
                    const int np = static_cast<int>(pop.size());
                    if (np < 3) continue;

                    // Sharing (uses current scaled fit values stored in fit[])
                    std::vector<double> sfit(np);
                    for (int i = 0; i < np; ++i) {
                        double sh = 0.0;
                        for (int j = 0; j < np; ++j) {
                            double d = euclidDist(pop[i], pop[j]);
                            if (d < m_sigma)
                                sh += 1.0 - std::pow(d / m_sigma, m_sharing_alpha);
                        }
                        sfit[i] = (sh > 1e-9) ? fit[i] / sh : fit[i];
                    }

                    std::uniform_int_distribution<int> pick(0, np - 1);
                    for (int i = 0; i < np; ++i) {
                        int r1, r2, r3;
                        do { r1 = pick(rng); } while (r1 == i);
                        do { r2 = pick(rng); } while (r2 == i || r2 == r1);
                        do { r3 = pick(rng); } while (r3 == i || r3 == r1 || r3 == r2);
                        std::vector<double> trial(D);
                        const int jr = static_cast<int>(u(rng) * D);
                        for (int d = 0; d < D; ++d) {
                            trial[d] = (u(rng) < PARTY_CR || d == jr)
                                ? clamp01(pop[r1][d] + m_de_F * (pop[r2][d] - pop[r3][d]))
                                : pop[i][d];
                        }
                        // accept/reject on scaled fitness
                        double tf = evalPartyScaled(trial, party);
                        if (tf > fit[i]) { pop[i] = trial; fit[i] = tf; }
                    }

                    // Clearing (using scaled fitness stored in fit[])
                    for (int i = 0; i < np; ++i) {
                        if (fit[i] <= DEAD + 1.0) continue;
                        for (int j = i + 1; j < np; ++j) {
                            if (fit[j] <= DEAD + 1.0) continue;
                            if (euclidDist(pop[i], pop[j]) < CLEARING_RADIUS) {
                                int loser = (sfit[i] >= sfit[j]) ? j : i;
                                for (auto& v : pop[loser]) v = u(rng);
                                fit[loser] = evalPartyScaled(pop[loser], party);
                            }
                        }
                    }
                }
            }

            // Extract per-party optima (store RAW fitness for logging/metrics)
            m_party_optima[party].clear();
            m_party_opt_fit[party].clear();
            for (size_t sp = 0; sp < m_competing_pops[party].size(); ++sp) {
                auto& fit = m_competing_fit[party][sp];
                auto  bfit = std::max_element(fit.begin(), fit.end());
                if (*bfit > DEAD + 1.0) {
                    size_t bi = static_cast<size_t>(std::distance(fit.begin(), bfit));
                    const auto& best_x = m_competing_pops[party][sp][bi];
                    m_party_optima[party].push_back(best_x);
                    m_party_opt_fit[party].push_back(evalParty(best_x, party)); // raw
                }
            }
        }

        //  Phase 2: CartesianSeed
        void CartesianSeed() {
            const int D = m_bench->dim();
            const int M = m_bench->numParties();
            const int npm = m_med_pop_size;
            std::mt19937 rng(static_cast<unsigned>(54321u));
            std::uniform_real_distribution<double> u(0.0, 1.0);

            std::vector<std::vector<double>> seeds;
            for (int p = 0; p < M; ++p)
                for (const auto& x : m_party_optima[p]) seeds.push_back(x);
            for (const auto& x0 : m_party_optima[0]) {
                for (const auto& x1 : m_party_optima[1]) {
                    std::vector<double> mid(D);
                    for (int d = 0; d < D; ++d) mid[d] = 0.5 * (x0[d] + x1[d]);
                    seeds.push_back(mid);
                }
            }

            m_mediating_pop.clear();
            m_mediating_fit.clear();
            m_mediating_obj.clear();

            for (const auto& s : seeds) {
                auto r = doEval(s);
                caInsert(s, r.consensus);
                osaInsert(s, r.consensus);
                if (static_cast<int>(m_mediating_pop.size()) >= npm) continue;
                bool close = false;
                for (const auto& ex : m_mediating_pop) {
                    if (euclidDist(s, ex) < m_sigma * 0.5) { close = true; break; }
                }
                if (!close) {
                    m_mediating_pop.push_back(s);
                    m_mediating_fit.push_back(r.consensus);
                    m_mediating_obj.push_back({ r.obj[0], r.obj[1] });
                }
            }

            while (static_cast<int>(m_mediating_pop.size()) < npm) {
                std::vector<double> x(D);
                for (auto& v : x) v = u(rng);
                auto r = doEval(x);
                caInsert(x, r.consensus);
                m_mediating_pop.push_back(x);
                m_mediating_fit.push_back(r.consensus);
                m_mediating_obj.push_back({ r.obj[0], r.obj[1] });
            }

            double best = *std::max_element(m_mediating_fit.begin(), m_mediating_fit.end());
            std::cout << "  Mediating pop seeded: " << m_mediating_pop.size()
                << " | best=" << std::fixed << std::setprecision(4) << best
                << " | ca=" << m_ca.size() << "\n";
        }

        //  Phase 3a: Evolve competing subpopulations
        void evolveCompeting() {
            const int M = m_bench->numParties();
            const int D = m_bench->dim();
            std::mt19937 rng(static_cast<unsigned>(m_generation * 12345u + 67890u));
            std::uniform_real_distribution<double> u(0.0, 1.0);

            for (int party = 0; party < M; ++party) {
                for (size_t sp = 0; sp < m_competing_pops[party].size(); ++sp) {
                    auto& pop = m_competing_pops[party][sp];
                    auto& fit = m_competing_fit[party][sp];
                    const int np = static_cast<int>(pop.size());
                    if (np < 3) continue;
                    std::uniform_int_distribution<int> pick(0, np - 1);
                    for (int i = 0; i < np; ++i) {
                        int r1, r2, r3;
                        do { r1 = pick(rng); } while (r1 == i);
                        do { r2 = pick(rng); } while (r2 == i || r2 == r1);
                        do { r3 = pick(rng); } while (r3 == i || r3 == r1 || r3 == r2);
                        std::vector<double> trial(D);
                        const int jr = static_cast<int>(u(rng) * D);
                        for (int d = 0; d < D; ++d) {
                            trial[d] = (u(rng) < m_de_CR || d == jr)
                                ? clamp01(pop[r1][d] + m_de_F * (pop[r2][d] - pop[r3][d]))
                                : pop[i][d];
                        }
                        // use scaled eval in competing loop too
                        double tf = evalPartyScaled(trial, party);
                        if (tf > fit[i]) { pop[i] = trial; fit[i] = tf; }
                    }

                    auto bfit_it = std::max_element(fit.begin(), fit.end());
                    if (*bfit_it > DEAD + 1.0) {
                        size_t bi = static_cast<size_t>(std::distance(fit.begin(), bfit_it));
                        auto r = doEval(pop[bi]); // raw eval for archive
                        caInsert(pop[bi], r.consensus);
                        osaInsert(pop[bi], r.consensus);
                        auto wt = std::min_element(m_mediating_fit.begin(), m_mediating_fit.end());
                        if (r.consensus > *wt) {
                            size_t wi = static_cast<size_t>(
                                std::distance(m_mediating_fit.begin(), wt));
                            m_mediating_pop[wi] = pop[bi];
                            m_mediating_fit[wi] = r.consensus;
                            m_mediating_obj[wi] = { r.obj[0], r.obj[1] };
                        }
                    }
                }

                m_party_optima[party].clear();
                m_party_opt_fit[party].clear();
                for (size_t sp = 0; sp < m_competing_pops[party].size(); ++sp) {
                    auto& fit = m_competing_fit[party][sp];
                    auto  bfit = std::max_element(fit.begin(), fit.end());
                    if (*bfit > DEAD + 1.0) {
                        size_t bi = static_cast<size_t>(std::distance(fit.begin(), bfit));
                        m_party_optima[party].push_back(m_competing_pops[party][sp][bi]);
                        m_party_opt_fit[party].push_back(
                            evalParty(m_competing_pops[party][sp][bi], party)); // raw
                    }
                }
            }
        }

        //  Phase 3b: Evolve mediating population
        //  C(x)-tournament DE + clearing
        void evolveMediating() {
            const int D = m_bench->dim();
            const int npm = static_cast<int>(m_mediating_pop.size());
            if (npm < 8) return;

            std::mt19937 rng(static_cast<unsigned>(m_generation * 444333u));
            std::uniform_real_distribution<double> u(0.0, 1.0);

            // Step 1: Inject per-party bests into weakest mediating slots
            for (int party = 0; party < m_bench->numParties(); ++party) {
                for (const auto& opt : m_party_optima[party]) {
                    auto r = doEval(opt);
                    caInsert(opt, r.consensus);
                    osaInsert(opt, r.consensus);
                    auto wt = std::min_element(m_mediating_fit.begin(), m_mediating_fit.end());
                    if (r.consensus > *wt) {
                        size_t wi = static_cast<size_t>(
                            std::distance(m_mediating_fit.begin(), wt));
                        m_mediating_pop[wi] = opt;
                        m_mediating_fit[wi] = r.consensus;
                        m_mediating_obj[wi] = { r.obj[0], r.obj[1] };
                    }
                }
            }

            // Step 2: C(x)-tournament DE/rand/1
            std::uniform_int_distribution<int> pick(0, npm - 1);
            for (int i = 0; i < npm; ++i) {
                int ta = pick(rng), tb = pick(rng);
                int r1 = (m_mediating_fit[ta] >= m_mediating_fit[tb]) ? ta : tb;
                int r2, r3;
                do { r2 = pick(rng); } while (r2 == r1);
                do { r3 = pick(rng); } while (r3 == r1 || r3 == r2);

                std::vector<double> trial(D);
                const int jr = static_cast<int>(u(rng) * D);
                for (int d = 0; d < D; ++d) {
                    trial[d] = (u(rng) < m_de_CR || d == jr)
                        ? clamp01(m_mediating_pop[r1][d]
                            + m_de_F * (m_mediating_pop[r2][d] - m_mediating_pop[r3][d]))
                        : m_mediating_pop[i][d];
                }

                auto tr = doEval(trial);
                if (tr.consensus > m_mediating_fit[i]) {
                    m_mediating_pop[i] = trial;
                    m_mediating_fit[i] = tr.consensus;
                    m_mediating_obj[i] = { tr.obj[0], tr.obj[1] };
                    caInsert(trial, tr.consensus);
                    osaInsert(trial, tr.consensus);
                }
            }

            // Step 3: Clearing pass
            for (int i = 0; i < npm; ++i) {
                if (m_mediating_fit[i] <= DEAD + 1.0) continue;
                for (int j = i + 1; j < npm; ++j) {
                    if (m_mediating_fit[j] <= DEAD + 1.0) continue;
                    if (euclidDist(m_mediating_pop[i], m_mediating_pop[j]) < CLEARING_RADIUS) {
                        int loser = (m_mediating_fit[i] >= m_mediating_fit[j]) ? j : i;
                        for (auto& v : m_mediating_pop[loser]) v = u(rng);
                        auto lr = doEval(m_mediating_pop[loser]);
                        m_mediating_fit[loser] = lr.consensus;
                        m_mediating_obj[loser] = { lr.obj[0], lr.obj[1] };
                        break;
                    }
                }
            }
        }

        //  Phase 3c: DM-oriented preservation
        void applyDMPreservation() {
            const int M = m_bench->numParties();
            const int npm = static_cast<int>(m_mediating_pop.size());
            const int D = m_bench->dim();

            for (int p = 0; p < M; ++p) {
                std::vector<std::pair<double, int>> scores;
                scores.reserve(npm);
                for (int i = 0; i < npm; ++i)
                    scores.push_back({ m_mediating_fit[i], i });
                std::sort(scores.begin(), scores.end(),
                    [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                        return a.first > b.first; });

                m_dm_elite[p].clear(); m_dm_elite_fit[p].clear();
                for (auto& [sc, idx] : scores) {
                    if (static_cast<int>(m_dm_elite[p].size()) >= m_K) break;
                    bool close = false;
                    for (const auto& ex : m_dm_elite[p]) {
                        if (euclidDist(m_mediating_pop[idx], ex) < m_sigma * 2.0) {
                            close = true; break;
                        }
                    }
                    if (!close) {
                        m_dm_elite[p].push_back(m_mediating_pop[idx]);
                        m_dm_elite_fit[p].push_back(sc);
                        osaInsert(m_mediating_pop[idx], sc);
                        caInsert(m_mediating_pop[idx], sc);
                    }
                }

                for (size_t ei = 0; ei < m_dm_elite[p].size(); ++ei) {
                    const auto& elite = m_dm_elite[p][ei];
                    bool present = false;
                    for (const auto& mp : m_mediating_pop) {
                        if (euclidDist(elite, mp) < 1e-6) { present = true; break; }
                    }
                    if (!present) {
                        auto r = doEval(elite);
                        auto wt = std::min_element(m_mediating_fit.begin(), m_mediating_fit.end());
                        if (r.consensus > *wt) {
                            size_t wi = static_cast<size_t>(
                                std::distance(m_mediating_fit.begin(), wt));
                            m_mediating_pop[wi] = elite;
                            m_mediating_fit[wi] = r.consensus;
                            m_mediating_obj[wi] = { r.obj[0], r.obj[1] };
                        }
                    }
                }
            }
        }

        //  Phase 3d: Stagnation detection + restart
        //
        //  On restart, also inject from uncovered basins.
        //  Uncovered basins (m_basin_ca_slot[b]=-1) were never explored
        //  force-inject a perturbed version of each into the mediating pop.
        bool checkAndRestart() {
            double cur = *std::max_element(
                m_mediating_fit.begin(), m_mediating_fit.end());
            m_conv_log.push_back({ m_generation, cur });

            if (cur > m_best_prev + 1e-6) { m_best_prev = cur; m_stagnation_cnt = 0; }
            else ++m_stagnation_cnt;
            if (m_stagnation_cnt < STAGNATION_LIMIT) return false;

            std::cout << "[RESTART gen=" << m_generation
                << " best=" << std::fixed << std::setprecision(4) << cur
                << " ca=" << m_ca.size() << "]\n";
            ++m_restart_count; m_stagnation_cnt = 0;

            const int D = m_bench->dim();
            const int M = m_bench->numParties();
            std::mt19937 rng(static_cast<unsigned>(
                m_generation * 55555u + (unsigned)m_restart_count * 9999u));
            std::normal_distribution<double>  N01(0.0, 1.0);
            std::uniform_real_distribution<double> u(0.0, 1.0);

            // Diversify competing pops
            for (int p = 0; p < M; ++p) {
                for (size_t sp = 0; sp < m_competing_pops[p].size(); ++sp) {
                    auto& pop = m_competing_pops[p][sp];
                    auto& fit = m_competing_fit[p][sp];
                    auto  bfit = std::max_element(fit.begin(), fit.end());
                    size_t bi = static_cast<size_t>(std::distance(fit.begin(), bfit));
                    for (size_t k = 0; k < pop.size(); ++k) {
                        if (k == bi) continue;
                        for (auto& v : pop[k]) v = u(rng);
                        fit[k] = evalPartyScaled(pop[k], p);
                    }
                }
            }

            // Re-seed mediating pop from Consensus Archive
            for (size_t i = 0; i < m_ca.size(); ++i) {
                auto wt = std::min_element(m_mediating_fit.begin(), m_mediating_fit.end());
                if (m_ca_fit[i] > *wt) {
                    size_t wi = static_cast<size_t>(
                        std::distance(m_mediating_fit.begin(), wt));
                    m_mediating_pop[wi] = m_ca[i];
                    m_mediating_fit[wi] = m_ca_fit[i];
                    auto r = doEval(m_ca[i]);
                    m_mediating_obj[wi] = { r.obj[0], r.obj[1] };
                }
            }

            // Force-inject from uncovered basins
            int injected_basins = 0;
            for (int b = 0; b < static_cast<int>(m_known_basins.size()); ++b) {
                if (m_basin_ca_slot[b] >= 0) continue; 
                // Perturb basin centre and inject into weakest mediating slot
                std::vector<double> seed = m_known_basins[b];
                const double perturb = 0.05;
                for (auto& v : seed) v = clamp01(v + perturb * N01(rng));
                auto r = doEval(seed);
                caInsert(seed, r.consensus);
                auto wt = std::min_element(m_mediating_fit.begin(), m_mediating_fit.end());
                size_t wi = static_cast<size_t>(
                    std::distance(m_mediating_fit.begin(), wt));
                m_mediating_pop[wi] = seed;
                m_mediating_fit[wi] = r.consensus;
                m_mediating_obj[wi] = { r.obj[0], r.obj[1] };
                ++injected_basins;
            }
            if (injected_basins > 0)
                std::cout << "  [Restart] Injected " << injected_basins
                << " uncovered-basin seeds\n";

            // Random exploration of underexplored regions
            const int n_replace = static_cast<int>(m_mediating_pop.size()) / 5;
            int inserted = 0;
            for (int att = 0; att < 5000 && inserted < n_replace; ++att) {
                std::vector<double> rsol(D);
                for (auto& v : rsol) v = u(rng);
                auto rr = doEval(rsol);
                if (rr.consensus < m_ca_min_quality) continue;
                bool too_close = false;
                for (const auto& ca_sol : m_ca) {
                    if (euclidDist(rsol, ca_sol) < 0.15) { too_close = true; break; }
                }
                if (!too_close) {
                    auto wt = std::min_element(m_mediating_fit.begin(), m_mediating_fit.end());
                    size_t wi = static_cast<size_t>(
                        std::distance(m_mediating_fit.begin(), wt));
                    m_mediating_pop[wi] = rsol;
                    m_mediating_fit[wi] = rr.consensus;
                    m_mediating_obj[wi] = { rr.obj[0], rr.obj[1] };
                    ++inserted;
                }
            }
            return true;
        }

        //  Phase 4: Joint refinement
        void JointRefine(std::vector<double>& x, unsigned seed = 0) {
            const int D = m_bench->dim();
            std::mt19937 rng(seed ? seed : static_cast<unsigned>(42u));
            std::normal_distribution<double> N01(0.0, 1.0);
            double cur_f = evalConsensus(x);
            for (int iter = 0; iter < m_refine_gens; ++iter) {
                double sig = std::max(0.005, 0.04 * std::exp(-3.0 * iter / m_refine_gens));
                std::vector<double> trial = x;
                for (int d = 0; d < D; ++d) trial[d] = clamp01(x[d] + sig * N01(rng));
                double tf = evalConsensus(trial);
                if (tf > cur_f) { x = trial; cur_f = tf; }
            }
        }

        //  Collect final candidates
        std::vector<std::vector<double>> collectCandidates() const {
            const int D = m_bench->dim();
            std::vector<std::vector<double>> cands = m_mediating_pop;
            for (int p = 0; p < m_bench->numParties(); ++p) {
                for (const auto& e : m_dm_elite[p])     cands.push_back(e);
                for (const auto& o : m_party_optima[p]) cands.push_back(o);
            }
            for (const auto& e : m_ca)  cands.push_back(e);
            for (const auto& e : m_osa) cands.push_back(e);

            std::vector<double> fits(cands.size());
            for (size_t i = 0; i < cands.size(); ++i)
                fits[i] = evalConsensus(cands[i]);

            std::vector<size_t> ord(cands.size());
            std::iota(ord.begin(), ord.end(), 0);
            std::sort(ord.begin(), ord.end(),
                [&](size_t a, size_t b) { return fits[a] > fits[b]; });

            std::vector<std::vector<double>> sel;
            for (size_t idx : ord) {
                bool close = false;
                for (const auto& ex : sel) {
                    if (euclidDist(cands[idx], ex) < 0.03) { close = true; break; }
                }
                if (!close) sel.push_back(cands[idx]);
            }
            return sel;
        }

        void adaptOsaRadius() {
            m_osa_radius = std::max(0.02, 0.04 * std::sqrt(2.0 / std::max(2, m_bench->dim())));
        }

        //  VISUALIZATION: landscape grids
        void saveLandscapes() const {
            const int D = m_bench->dim();
            const int RES = VIS_RES;
            std::filesystem::create_directories(m_out_dir);

            struct Sample { double x0, x1, f0, f1, cons; };
            std::vector<Sample> grid;
            grid.reserve(static_cast<size_t>(RES * RES));

            for (int i = 0; i < RES; ++i) {
                for (int j = 0; j < RES; ++j) {
                    std::vector<double> x(D, 0.5);
                    x[0] = static_cast<double>(i) / (RES - 1);
                    if (D >= 2) x[1] = static_cast<double>(j) / (RES - 1);
                    auto r = doEval(x);
                    grid.push_back({ x[0], (D >= 2 ? x[1] : 0.5), r.obj[0], r.obj[1], r.consensus });
                }
            }

            std::vector<int> prank(grid.size(), 0);
            for (size_t i = 0; i < grid.size(); ++i) {
                for (size_t j = 0; j < grid.size(); ++j) {
                    if (i == j) continue;
                    bool j_dom_i = (grid[j].f0 >= grid[i].f0 - 1e-9) &&
                        (grid[j].f1 >= grid[i].f1 - 1e-9) &&
                        (grid[j].f0 > grid[i].f0 + 1e-9 ||
                            grid[j].f1 > grid[i].f1 + 1e-9);
                    if (j_dom_i) ++prank[i];
                }
            }

            std::string dim_note = " D=" + std::to_string(D)
                + (D > 2 ? " (slice: x0-x1, others=0.5)" : "");

            auto write = [&](const std::string& fname,
                const std::string& header,
                std::function<double(size_t)> val) {
                    std::ofstream f(m_out_dir + "/" + fname);
                    if (!f.is_open()) { std::cerr << "[WARN] Cannot write " << fname << "\n"; return; }
                    f << header << dim_note << "\n# x0 x1 value\n"
                        << std::fixed << std::setprecision(6);
                    for (size_t k = 0; k < grid.size(); ++k)
                        f << grid[k].x0 << " " << grid[k].x1 << " " << val(k) << "\n";
                };

            write("landscape_f0.txt", "# Party-0 landscape F0(x)", [&](size_t k) { return grid[k].f0; });
            write("landscape_f1.txt", "# Party-1 landscape F1(x)", [&](size_t k) { return grid[k].f1; });
            write("landscape_consensus.txt", "# Consensus C(x)=min(F0,F1)", [&](size_t k) { return grid[k].cons; });
            write("landscape_rank.txt", "# Two-objective rank landscape", [&](size_t k) { return static_cast<double>(prank[k]); });

            std::cout << "  [VIS] Landscape grids saved (" << RES << "x" << RES
                << ", D=" << D << ")\n";
        }

        void saveOptima() const {
            std::filesystem::create_directories(m_out_dir);

            {
                std::ofstream f(m_out_dir + "/shared_optima.txt");
                if (f.is_open()) {
                    f << "# Ground-truth shared optima\n# x0..x(D-1) f0 f1 consensus\n"
                        << std::fixed << std::setprecision(6);
                    for (const auto& o : m_bench->sharedOptima()) {
                        std::vector<double> xd(o.begin(), o.end());
                        auto r = doEval(xd);
                        for (double v : xd) f << v << " ";
                        f << r.obj[0] << " " << r.obj[1] << " " << r.consensus << "\n";
                    }
                }
            }

            for (int p = 0; p < m_bench->numParties(); ++p) {
                std::ofstream f(m_out_dir + "/party" + std::to_string(p) + "_optima.txt");
                if (!f.is_open()) continue;
                f << "# Party-" << p << " optima (Phase 1b)\n"
                    << "# x0..x(D-1) party_fit f0 f1 consensus\n"
                    << std::fixed << std::setprecision(6);
                for (size_t k = 0; k < m_party_optima[p].size(); ++k) {
                    auto r = doEval(m_party_optima[p][k]);
                    for (double v : m_party_optima[p][k]) f << v << " ";
                    f << m_party_opt_fit[p][k] << " "
                        << r.obj[0] << " " << r.obj[1] << " " << r.consensus << "\n";
                }
            }
        }

        void saveVisualizationData(int tag) const {
            std::filesystem::create_directories(m_out_dir);
            const int M = m_bench->numParties();

            {
                std::ofstream f(m_out_dir + "/gen_" + std::to_string(tag) + "_mediating.txt");
                if (f.is_open()) {
                    f << "# gen=" << tag << " mediating population\n# x0..x(D-1) f0 f1 consensus\n"
                        << std::fixed << std::setprecision(6);
                    for (size_t i = 0; i < m_mediating_pop.size(); ++i) {
                        if (m_mediating_fit[i] <= DEAD + 1.0) continue;
                        for (double v : m_mediating_pop[i]) f << v << " ";
                        if (i < m_mediating_obj.size())
                            f << m_mediating_obj[i][0] << " " << m_mediating_obj[i][1] << " ";
                        f << m_mediating_fit[i] << "\n";
                    }
                }
            }

            for (int p = 0; p < M; ++p) {
                std::ofstream f(m_out_dir + "/gen_" + std::to_string(tag)
                    + "_party" + std::to_string(p) + ".txt");
                if (!f.is_open()) continue;
                f << "# gen=" << tag << " party=" << p << "\n# sp x0..x(D-1) party_fit\n"
                    << std::fixed << std::setprecision(6);
                for (size_t sp = 0; sp < m_competing_pops[p].size(); ++sp)
                    for (size_t i = 0; i < m_competing_pops[p][sp].size(); ++i) {
                        if (m_competing_fit[p][sp][i] <= DEAD + 1.0) continue;
                        f << sp;
                        for (double v : m_competing_pops[p][sp][i]) f << " " << v;
                        f << " " << m_competing_fit[p][sp][i] << "\n";
                    }
            }

            {
                std::ofstream f(m_out_dir + "/gen_" + std::to_string(tag) + "_archive.txt");
                if (f.is_open()) {
                    f << "# gen=" << tag << " CA size=" << m_ca.size() << "\n# x0..x(D-1) consensus\n"
                        << std::fixed << std::setprecision(6);
                    for (size_t i = 0; i < m_ca.size(); ++i) {
                        for (double v : m_ca[i]) f << v << " ";
                        f << m_ca_fit[i] << "\n";
                    }
                }
            }
        }

    public:
        //  Constructor
        //  m_ca_min_quality is dimension-dependent.
        MPM_CoEA(std::shared_ptr<MPMMO_Benchmark> bench,
            std::shared_ptr<Environment>     env,
            int    K = 8,
            int    niche_pop = 20,
            int    med_pop = 300,
            int    party_gens = 120,
            int    coevo_gens = 400,
            int    refine_gens = 80,
            double sigma = 0.12,
            double de_F = 0.5,
            double de_CR = 0.9,
            double sharing_a = 2.0,
            const std::string& out_dir = "Visualization/solutions")
            : m_bench(bench), m_env_keeper(env),
            m_K(K), m_niche_pop(niche_pop), m_med_pop_size(med_pop),
            m_party_gens(party_gens), m_coevo_gens(coevo_gens),
            m_refine_gens(refine_gens), m_sigma(sigma),
            m_de_F(de_F), m_de_CR(de_CR), m_sharing_alpha(sharing_a),
            m_osa_radius(0.04), m_out_dir(out_dir)
        {
            // dim-dependent CA quality threshold
            const int D = bench->dim();
            if (D <= 2)       m_ca_min_quality = 20.0;
            else if (D <= 5)  m_ca_min_quality = 5.0;
            else              m_ca_min_quality = 1.0;

            std::filesystem::create_directories(out_dir);
            const int M = bench->numParties();
            m_party_optima.resize(M);
            m_party_opt_fit.resize(M);
            m_competing_pops.resize(M);
            m_competing_fit.resize(M);
            m_dm_elite.resize(M);
            m_dm_elite_fit.resize(M);
        }

        //  Main run
        void run() {
            const int M = m_bench->numParties();
            const int D = m_bench->dim();

            std::cout << "  [VIS] Generating landscape grids...\n";
            saveLandscapes();

            // Phase 1a: NBC
            std::cout << "\n>>> [MPM-CoEA] Phase 1a: NBC init\n";
            for (int p = 0; p < M; ++p) NBC_init(p);
            for (int p = 0; p < M; ++p)
                std::cout << "  P" << p << ": "
                << m_party_optima[p].size() << " initial optima\n";

            // Phase 1b: ClearingDE
            std::cout << ">>> [MPM-CoEA] Phase 1b: ClearingDE ("
                << m_party_gens << " gens, scaled eval D=" << D << ")\n";
            for (int p = 0; p < M; ++p) ClearingDE(p);
            for (int p = 0; p < M; ++p)
                for (size_t k = 0; k < m_party_optima[p].size(); ++k)
                    std::cout << "  [P" << p << " opt" << k << "] f="
                    << std::fixed << std::setprecision(4)
                    << m_party_opt_fit[p][k] << "\n";

            // Build basin registry from Phase 1b results
            buildKnownBasins();

            saveOptima();

            // Phase 2
            std::cout << ">>> [MPM-CoEA] Phase 2: CartesianSeed\n";
            CartesianSeed();

            // Phase 3
            std::cout << ">>> [MPM-CoEA] Phase 3: Co-evolution ("
                << m_coevo_gens << " gens)\n";
            std::cout << "  CA_MIN_QUALITY=" << m_ca_min_quality
                << "  Basins=" << m_known_basins.size() << "\n";

            m_best_prev = *std::max_element(
                m_mediating_fit.begin(), m_mediating_fit.end());
            m_conv_log.push_back({ 0, m_best_prev });
            std::cout << "  Gen    0 | best="
                << std::fixed << std::setprecision(4) << m_best_prev
                << " | stag=0 | ca=" << m_ca.size() << "\n";
            saveVisualizationData(0);

            static const int snaps[] = { 10, 50, 100, 200, 300, 400, -1 };

            for (m_generation = 1; m_generation <= m_coevo_gens; ++m_generation) {
                adaptOsaRadius();
                evolveCompeting();
                evolveMediating();
                applyDMPreservation();
                checkAndRestart();

                for (int k = 0; snaps[k] > 0; ++k)
                    if (m_generation == snaps[k] || m_generation == m_coevo_gens)
                        saveVisualizationData(m_generation);

                if (m_generation % 5 == 0) {
                    double best = *std::max_element(
                        m_mediating_fit.begin(), m_mediating_fit.end());
                    // Count covered basins
                    int covered = 0;
                    for (int b : m_basin_ca_slot) if (b >= 0) ++covered;
                    std::cout << "  Gen " << std::setw(4) << m_generation
                        << " | best=" << std::fixed << std::setprecision(4) << best
                        << " | stag=" << m_stagnation_cnt
                        << " | rst=" << m_restart_count
                        << " | ca=" << m_ca.size()
                        << " | cov=" << covered << "/" << m_known_basins.size()
                        << " | osa=" << m_osa.size() << "\n";
                }
            }

            // Phase 4: Refinement
            std::cout << ">>> [MPM-CoEA] Phase 4: JointRefine\n";
            auto finals = collectCandidates();
            finals.resize(std::min(static_cast<int>(finals.size()), 120));
            for (size_t i = 0; i < finals.size(); ++i) {
                JointRefine(finals[i], static_cast<unsigned>(i * 13u + 42u));
                double rf = evalConsensus(finals[i]);
                caInsert(finals[i], rf);
                osaInsert(finals[i], rf);
            }
            for (size_t i = 0; i < finals.size() && i < m_mediating_pop.size(); ++i) {
                auto r = doEval(finals[i]);
                m_mediating_pop[i] = finals[i];
                m_mediating_fit[i] = r.consensus;
                m_mediating_obj[i] = { r.obj[0], r.obj[1] };
            }

            int covered = 0;
            for (int b : m_basin_ca_slot) if (b >= 0) ++covered;
            double best_f = *std::max_element(
                m_mediating_fit.begin(), m_mediating_fit.end());
            std::cout << "  Best (post-refine): "
                << std::fixed << std::setprecision(4) << best_f
                << "  CA=" << m_ca.size()
                << "  BasinsCovered=" << covered << "/" << m_known_basins.size()
                << "  OSA=" << m_osa.size() << "\n";
        }

        //  Accessors
        double getBestFitness() const {
            if (m_mediating_fit.empty()) return DEAD;
            return *std::max_element(m_mediating_fit.begin(), m_mediating_fit.end());
        }
        std::vector<double> getBestSolution() const {
            if (m_mediating_pop.empty()) return {};
            auto it = std::max_element(m_mediating_fit.begin(), m_mediating_fit.end());
            return m_mediating_pop[
                static_cast<size_t>(std::distance(m_mediating_fit.begin(), it))];
        }
        std::vector<std::vector<double>> getAllCandidates() const {
            return collectCandidates();
        }
        const std::vector<std::vector<double>>& getConsensusArchive()    const { return m_ca; }
        const std::vector<double>& getConsensusArchiveFit() const { return m_ca_fit; }

        //  Save functions
        void saveConvergenceCurve() const {
            std::filesystem::create_directories(m_out_dir);
            std::ofstream out(m_out_dir + "/convergence.txt");
            if (!out.is_open()) return;
            out << "# generation best_consensus ca_size\n";
            for (const auto& p : m_conv_log)
                out << p.first << " "
                << std::fixed << std::setprecision(6) << p.second
                << " " << m_ca.size() << "\n";
        }

        void saveConsensusArchive() const {
            std::filesystem::create_directories(m_out_dir);
            std::ofstream out(m_out_dir + "/consensus_archive.txt");
            if (!out.is_open()) return;
            out << "# Consensus Archive (basin-tagged)\n"
                << "# x0..x(D-1) f0 f1 consensus\n"
                << std::fixed << std::setprecision(6);
            for (size_t i = 0; i < m_ca.size(); ++i) {
                for (double v : m_ca[i]) out << v << " ";
                auto r = doEval(m_ca[i]);
                out << r.obj[0] << " " << r.obj[1] << " " << m_ca_fit[i] << "\n";
            }
        }

        void saveFinalMediatingPop() const {
            std::filesystem::create_directories(m_out_dir);
            std::ofstream out(m_out_dir + "/final_mediating_pop.txt");
            if (!out.is_open()) return;
            out << "# Final mediating pop (post Phase 4)\n# x0..x(D-1) f0 f1 consensus\n"
                << std::fixed << std::setprecision(6);
            for (size_t i = 0; i < m_mediating_pop.size(); ++i) {
                if (m_mediating_fit[i] <= DEAD + 1.0) continue;
                for (double v : m_mediating_pop[i]) out << v << " ";
                if (i < m_mediating_obj.size())
                    out << m_mediating_obj[i][0] << " " << m_mediating_obj[i][1] << " ";
                out << m_mediating_fit[i] << "\n";
            }
        }
    };  // class MPM_CoEA

    //  Benchmark result struct
    struct BenchmarkResult {
        std::string name;
        int         suite_id = 0;
        int         dim = 0;
        int         nopt = 0;
        double      sr_strict = 0, anof_strict = 0, mpr_strict = 0;
        double      sr_loose = 0, anof_loose = 0, mpr_loose = 0;
        int         ca_size = 0;
        double      best_fitness = 0;
    };

    //  Single problem instance runner
    inline BenchmarkResult runSingleBenchmark(
        int         suite_id,
        size_t      dimension,
        const std::string& vis_dir,
        int         coevo_gens = 400)
    {
        std::cout << "\n[runSingleBenchmark] P" << suite_id
            << "  D=" << dimension << "\n";

        std::shared_ptr<Environment> env(
            generateEnvironmentByFactory("free_peaks_multiparty"));
        env->setProblem(generateProblemByFactory("free_peaks_multiparty"));

        ParameterMap pm;
        pm["suite_id"] = suite_id;
        pm["problem_dimension"] = dimension;
        env->problem()->inputParameters().input(pm);
        env->initializeProblem(0.5);

        auto* fp = dynamic_cast<FreePeaksMultiParty*>(env->problem());
        if (!fp) throw std::runtime_error("FreePeaksMultiParty cast failed");

        const auto& spec = fp->currentSpec();
        const int nopt = static_cast<int>(spec.shared_optima.size());
        std::cout << "  problem=" << spec.name
            << "  dim=" << dimension
            << "  nopt=" << nopt
            << "  feature=" << spec.feature << "\n";

        auto bench = std::make_shared<MPMMO_Benchmark>(fp, env.get());

        MPM_CoEA solver(bench, env,
            /*K*/           8,
            /*niche_pop*/  20,
            /*med_pop*/   300,
            /*party_gens*/120,
            /*coevo_gens*/coevo_gens,
            /*refine_gens*/ 80,
            /*sigma*/      0.12,
            /*de_F*/       0.5,
            /*de_CR*/      0.9,
            /*sharing_a*/  2.0,
            vis_dir);

        solver.run();
        solver.saveConvergenceCurve();
        solver.saveConsensusArchive();
        solver.saveFinalMediatingPop();

        auto candidates = solver.getAllCandidates();
        auto ms = evaluateCEC2015Metrics(candidates, *bench, 0.01, 1.0);
        auto ml = evaluateCEC2015Metrics(candidates, *bench, 0.05, 2.0);

        std::cout << "  CEC-2015 strict (eps_x=0.01, eps_f=1.0):\n"
            << "    SR=" << ms.SR << "%  ANOF=" << ms.ANOF
            << "  MPR=" << ms.MPR << "\n"
            << "  CEC-2015 loose (eps_x=0.05, eps_f=2.0):\n"
            << "    SR=" << ml.SR << "%  ANOF=" << ml.ANOF
            << "  MPR=" << ml.MPR << "\n"
            << "  CA=" << solver.getConsensusArchive().size()
            << "  best_C=" << std::fixed << std::setprecision(4)
            << solver.getBestFitness() << "\n";

        {
            std::filesystem::create_directories(vis_dir);
            std::ofstream mf(vis_dir + "/cec2015_metrics.txt");
            if (mf.is_open()) {
                mf << std::fixed << std::setprecision(4)
                    << "# CEC-2015 metrics  P" << suite_id << "  D=" << dimension << "\n"
                    << "suite_id       " << suite_id << "\n"
                    << "problem_name   " << spec.name << "\n"
                    << "feature        " << spec.feature << "\n"
                    << "dim            " << dimension << "\n"
                    << "nopt           " << nopt << "\n"
                    << "# strict eps_x=0.01 eps_f=1.0\n"
                    << "SR_strict      " << ms.SR << "\n"
                    << "ANOF_strict    " << ms.ANOF << "\n"
                    << "MPR_strict     " << ms.MPR << "\n"
                    << "# loose  eps_x=0.05 eps_f=2.0\n"
                    << "SR_loose       " << ml.SR << "\n"
                    << "ANOF_loose     " << ml.ANOF << "\n"
                    << "MPR_loose      " << ml.MPR << "\n"
                    << "ca_size        " << solver.getConsensusArchive().size() << "\n"
                    << "best_consensus " << solver.getBestFitness() << "\n";
            }
        }

        BenchmarkResult r;
        r.name = spec.name;
        r.suite_id = suite_id;
        r.dim = static_cast<int>(dimension);
        r.nopt = nopt;
        r.sr_strict = ms.SR;    r.anof_strict = ms.ANOF;  r.mpr_strict = ms.MPR;
        r.sr_loose = ml.SR;    r.anof_loose = ml.ANOF;  r.mpr_loose = ml.MPR;
        r.ca_size = static_cast<int>(solver.getConsensusArchive().size());
        r.best_fitness = solver.getBestFitness();
        return r;
    }

    //  Full suite runner
    //  Adaptive generation according to dimension.
    inline void runBenchmarkSuiteAll() {
        registerInstance();

        const std::string vis_base =
            "E:\\HITSZ\\Research\\Multimodal_Multiparty_Optimization"
            "\\ThesisProject\\Visualization";

        std::cout << " MPM-CoEA Benchmark Suite\n";

        const std::vector<int>    suite_ids = { 1,2,3,4,5,6,7,8,9,10,11,12 };
        const std::vector<size_t> dimensions = { 2, 3, 5, 10 };

        struct Inst { int sid; size_t dim; int gens; };
        std::vector<Inst> instances;
        for (int sid : suite_ids) {
            for (size_t d : dimensions) {
                // Adaptive budget - 2D converges fast, high-D needs more time
                int gens;
                if (d <= 2)  gens = 200;
                else if (d <= 5)  gens = 500;
                else              gens = 800;
                instances.push_back({ sid, d, gens });
            }
        }

        std::vector<BenchmarkResult> results;
        int idx = 0;
        for (const auto& inst : instances) {
            ++idx;
            std::cout << "\n[" << idx << "/" << instances.size() << "]  "
                << "P" << inst.sid << "  D=" << inst.dim
                << "  gens=" << inst.gens << "\n";
            std::string vis_dir = vis_base + "/mpmcoea_P"
                + std::to_string(inst.sid) + "_D" + std::to_string(inst.dim);
            try {
                results.push_back(runSingleBenchmark(inst.sid, inst.dim, vis_dir, inst.gens));
            }
            catch (const std::exception& ex) {
                std::cerr << "  [FAIL] P" << inst.sid
                    << " D=" << inst.dim << ": " << ex.what() << "\n";
            }
        }

        // Summary table
        std::cout << "\n>>> Results (" << results.size()
            << "/" << instances.size() << ")\n";
        std::cout << std::left
            << std::setw(32) << "Problem"
            << std::setw(5) << "Sid"
            << std::setw(5) << "Dim"
            << std::setw(5) << "Opt"
            << std::setw(9) << "SR(s)%"
            << std::setw(9) << "ANOF(s)"
            << std::setw(9) << "SR(l)%"
            << std::setw(9) << "ANOF(l)"
            << std::setw(7) << "CA_sz"
            << std::setw(10) << "Best_C"
            << "\n" << std::string(100, '-') << "\n";
        for (const auto& r : results) {
            std::cout << std::left
                << std::setw(32) << r.name
                << std::setw(5) << r.suite_id
                << std::setw(5) << r.dim
                << std::setw(5) << r.nopt
                << std::setw(9) << std::fixed << std::setprecision(1) << r.sr_strict
                << std::setw(9) << std::setprecision(3) << r.anof_strict
                << std::setw(9) << std::setprecision(1) << r.sr_loose
                << std::setw(9) << std::setprecision(3) << r.anof_loose
                << std::setw(7) << r.ca_size
                << std::setw(10) << std::setprecision(2) << r.best_fitness
                << "\n";
        }

        // CSV
        const std::string csv_path = vis_base + "/mpmcoea_suite_summary.csv";
        {
            std::ofstream csv(csv_path);
            if (csv.is_open()) {
                csv << "problem_name,suite_id,dim,nopt,"
                    "sr_strict,anof_strict,mpr_strict,"
                    "sr_loose,anof_loose,mpr_loose,"
                    "ca_size,best_fitness\n"
                    << std::fixed << std::setprecision(4);
                for (const auto& r : results)
                    csv << r.name << "," << r.suite_id << ","
                    << r.dim << "," << r.nopt << ","
                    << r.sr_strict << "," << r.anof_strict << "," << r.mpr_strict << ","
                    << r.sr_loose << "," << r.anof_loose << "," << r.mpr_loose << ","
                    << r.ca_size << "," << r.best_fitness << "\n";
                std::cout << "[INFO] CSV: " << csv_path << "\n";
            }
        }
        std::cout << "\n>>> [MPM-CoEA Suite] Complete.\n";
    }

    //  Entry point
    inline void runMPMCoEAExperiment() {
        runBenchmarkSuiteAll();
    }

}  // namespace ofec

#endif  // OFEC_MPMCOEA_SOLVER_HPP