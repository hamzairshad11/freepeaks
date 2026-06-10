#ifndef OFEC_MULTIPARTY_NBN_CMAES_SOLVER_HPP
#define OFEC_MULTIPARTY_NBN_CMAES_SOLVER_HPP

#include "interface.h"
#include "../core/problem/solution.h"
#include "../instance/problem/continuous/free_peaks/free_peaks_multiparty.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

namespace ofec {

namespace multiparty_nbn_cmaes {

struct Individual {
    std::vector<Real> x;
    std::array<Real, 2> obj { 0, 0 };
    Real consensus = 0;
};

struct NbdEdge {
    int from = -1;
    int to = -1;
    Real length = 0;
};

struct Cluster {
    std::vector<int> members;
    int best = -1;
};

struct Config {
    int iterations = 220;
    int party_pop_size = 180;
    int coord_pop_size = 220;
    int max_clusters = 12;
    int niche_lambda = 18;
    int coord_children = 120;
    int consensus_archive_max = 80;
    int nbn_source_cap = 1800;
    int max_evaluations = 1000000;
    Real archive_radius = 0.035;
    Real nbn_source_min_dist = 0.006;
};

inline Real sqr(Real v) { return v * v; }

inline Real dist2(const std::vector<Real>& a, const std::vector<Real>& b) {
    Real value = 0;
    const size_t n = std::min(a.size(), b.size());
    for (size_t i = 0; i < n; ++i) value += sqr(a[i] - b[i]);
    return value;
}

inline Real dist(const std::vector<Real>& a, const std::vector<Real>& b) {
    return std::sqrt(dist2(a, b));
}

inline Real clamp01(Real v) {
    return std::max<Real>(0, std::min<Real>(1, v));
}

class Solver {
public:
    Solver(FreePeaksMultiParty* problem, Environment* env, int suite_id,
        const std::string& output_dir, const Config& config = Config())
        : m_problem(problem), m_env(env), m_suite_id(suite_id), m_output_dir(output_dir),
          m_dim(problem ? problem->numberVariables() : 2),
          m_rng(static_cast<unsigned>(20250609u + suite_id * 131u + m_dim * 17u)),
          m_config(config) {}

    void run() {
        std::filesystem::create_directories(m_output_dir);
        initializePopulations();

        std::ofstream log(m_output_dir + "/convergence.txt");
        log << "# iter evaluations best_consensus history_archive_size consensus_archive_size party0_clusters party1_clusters\n";

        for (int iter = 0; iter < m_config.iterations && remainingEvaluations() > 0; ++iter) {
            updateArchives();
            std::array<std::vector<Individual>, 2> nbn_sources;
            std::array<std::vector<Cluster>, 2> clusters;
            for (int party = 0; party < 2; ++party) {
                if (remainingEvaluations() <= 0) break;
                nbn_sources[party] = selectNbnSource(party);
                clusters[party] = nbdClusters(nbn_sources[party], party);
                rebuildPartyPopFromClusters(party, clusters[party], nbn_sources[party]);
                evolvePartyByClusters(party, clusters[party], nbn_sources[party]);
            }

            updateArchives();
            if (remainingEvaluations() > 0) {
                updateCoordinator(clusters, nbn_sources);
                evolveCoordinator();
            }
            updateArchives();

            if (iter % 10 == 0 || iter + 1 == m_config.iterations || remainingEvaluations() <= 0) {
                log << iter << ' ' << m_evaluations << ' ' << std::fixed << std::setprecision(10)
                    << bestConsensus() << ' ' << m_archive.size() << ' '
                    << m_consensus_archive.size() << ' '
                    << clusters[0].size() << ' ' << clusters[1].size() << "\n";
                std::cout << "[suite " << m_suite_id << "] iter " << iter
                          << " best=" << std::fixed << std::setprecision(4)
                          << bestConsensus()
                          << " FE=" << m_evaluations
                          << " history=" << m_archive.size()
                          << " consensus_archive=" << m_consensus_archive.size()
                          << " clusters=(" << clusters[0].size() << "," << clusters[1].size() << ")\n";
            }
        }

        savePopulation(m_output_dir + "/party0_population.txt", m_party_pop[0]);
        savePopulation(m_output_dir + "/party1_population.txt", m_party_pop[1]);
        savePopulation(m_output_dir + "/coordination_population.txt", m_coord_pop);
        savePopulation(m_output_dir + "/archive_all_history.txt", m_archive);
        savePopulation(m_output_dir + "/archive_consensus_optima.txt", m_consensus_archive);
        saveNbdEdges(m_output_dir + "/party0_nbn_edges.txt", selectNbnSource(0), 0);
        saveNbdEdges(m_output_dir + "/party1_nbn_edges.txt", selectNbnSource(1), 1);
        saveRunInfo(m_output_dir + "/run_info.txt");
    }

private:
    int remainingEvaluations() const {
        return std::max(0, m_config.max_evaluations - m_evaluations);
    }

    Individual evaluate(const std::vector<Real>& x) {
        auto sol = std::unique_ptr<SolutionBase>(m_problem->createSolution());
        auto& var = dynamic_cast<Solution<VariableVector<Real>>&>(*sol).variable().vector();
        for (size_t d = 0; d < m_dim; ++d) var[d] = clamp01(x[d]);
        sol->evaluate(m_env, false);
        Individual ind;
        ind.x.assign(var.begin(), var.end());
        ind.obj = { sol->objective(0), sol->objective(1) };
        ind.consensus = std::min(ind.obj[0], ind.obj[1]);
        ++m_evaluations;
        return ind;
    }

    void initializePopulations() {
        std::uniform_real_distribution<Real> U(0, 1);
        for (int party = 0; party < 2; ++party) {
            m_party_pop[party].clear();
            for (int i = 0; i < m_config.party_pop_size && remainingEvaluations() > 0; ++i) {
                std::vector<Real> x(m_dim);
                for (auto& v : x) v = U(m_rng);
                m_party_pop[party].push_back(evaluate(x));
            }
        }
        m_coord_pop.clear();
        for (int i = 0; i < m_config.coord_pop_size && remainingEvaluations() > 0; ++i) {
            std::vector<Real> x(m_dim);
            for (auto& v : x) v = U(m_rng);
            m_coord_pop.push_back(evaluate(x));
        }
        updateArchives();
    }

    std::vector<Individual> selectNbnSource(int party) const {
        std::vector<Individual> source = m_archive;
        std::sort(source.begin(), source.end(), [&](const Individual& a, const Individual& b) {
            if (std::abs(a.obj[party] - b.obj[party]) > 1e-12) return a.obj[party] > b.obj[party];
            return a.consensus > b.consensus;
        });
        diverseTruncate(source, m_config.nbn_source_cap, m_config.nbn_source_min_dist);
        return source;
    }

    std::vector<NbdEdge> nearestBetterEdges(const std::vector<Individual>& pop, int party) const {
        std::vector<NbdEdge> edges;
        edges.reserve(pop.size());
        for (int i = 0; i < static_cast<int>(pop.size()); ++i) {
            int best_j = -1;
            Real best_d = std::numeric_limits<Real>::max();
            for (int j = 0; j < static_cast<int>(pop.size()); ++j) {
                if (pop[j].obj[party] <= pop[i].obj[party] + 1e-12) continue;
                Real d = dist(pop[i].x, pop[j].x);
                if (d < best_d) {
                    best_d = d;
                    best_j = j;
                }
            }
            if (best_j >= 0) edges.push_back({ i, best_j, best_d });
        }
        return edges;
    }

    std::vector<Cluster> nbdClusters(const std::vector<Individual>& pop, int party) const {
        const int n = static_cast<int>(pop.size());
        auto edges = nearestBetterEdges(pop, party);
        if (n == 0) return {};
        if (edges.empty()) {
            std::vector<Cluster> clusters;
            for (int i = 0; i < n; ++i) clusters.push_back({ { i }, i });
            return clusters;
        }

        Real mean = 0;
        for (const auto& e : edges) mean += e.length;
        mean /= std::max<int>(1, static_cast<int>(edges.size()));
        Real var = 0;
        for (const auto& e : edges) var += sqr(e.length - mean);
        var /= std::max<int>(1, static_cast<int>(edges.size()));
        const Real threshold = mean + 0.8 * std::sqrt(var);

        std::vector<std::vector<int>> adj(n);
        for (const auto& e : edges) {
            if (e.length <= threshold) {
                adj[e.from].push_back(e.to);
                adj[e.to].push_back(e.from);
            }
        }

        std::vector<int> seen(n, 0);
        std::vector<Cluster> clusters;
        for (int i = 0; i < n; ++i) {
            if (seen[i]) continue;
            Cluster c;
            std::queue<int> q;
            q.push(i);
            seen[i] = 1;
            while (!q.empty()) {
                int u = q.front(); q.pop();
                c.members.push_back(u);
                for (int v : adj[u]) if (!seen[v]) {
                    seen[v] = 1;
                    q.push(v);
                }
            }
            c.best = *std::max_element(c.members.begin(), c.members.end(),
                [&](int a, int b) { return pop[a].obj[party] < pop[b].obj[party]; });
            clusters.push_back(std::move(c));
        }
        std::sort(clusters.begin(), clusters.end(), [&](const Cluster& a, const Cluster& b) {
            return pop[a.best].obj[party] > pop[b.best].obj[party];
        });
        if (static_cast<int>(clusters.size()) > m_config.max_clusters)
            clusters.resize(m_config.max_clusters);
        return clusters;
    }

    void rebuildPartyPopFromClusters(int party, const std::vector<Cluster>& clusters,
        const std::vector<Individual>& source) {
        std::vector<Individual> rebuilt;
        rebuilt.reserve(m_config.party_pop_size + static_cast<int>(clusters.size()));

        for (const auto& c : clusters) {
            if (c.best >= 0 && c.best < static_cast<int>(source.size()))
                rebuilt.push_back(source[c.best]);
            std::vector<int> members = c.members;
            std::sort(members.begin(), members.end(), [&](int a, int b) {
                return source[a].obj[party] > source[b].obj[party];
            });
            for (int idx : members) {
                if (idx >= 0 && idx < static_cast<int>(source.size()))
                    rebuilt.push_back(source[idx]);
                if (static_cast<int>(rebuilt.size()) >= m_config.party_pop_size * 2) break;
            }
        }

        if (rebuilt.empty()) rebuilt = source;
        truncateByFitnessAndDiversity(rebuilt, party, m_config.party_pop_size, 0.012);
        m_party_pop[party].swap(rebuilt);
    }

    void evolvePartyByClusters(int party, const std::vector<Cluster>& clusters,
        const std::vector<Individual>& source) {
        std::normal_distribution<Real> N(0, 1);
        std::vector<Individual> offspring;
        offspring.reserve(static_cast<size_t>(clusters.size() * m_config.niche_lambda));

        for (const auto& c : clusters) {
            if (c.best < 0) continue;
            const auto center = source[c.best].x;
            Real spread = 0.02;
            for (int idx : c.members) spread += dist(center, source[idx].x);
            spread /= std::max<int>(1, static_cast<int>(c.members.size()));
            Real sigma = std::max<Real>(0.015, std::min<Real>(0.18, 0.65 * spread + 0.015));

            for (int k = 0; k < m_config.niche_lambda && remainingEvaluations() > 0; ++k) {
                std::vector<Real> y(m_dim);
                for (size_t d = 0; d < m_dim; ++d) y[d] = clamp01(center[d] + sigma * N(m_rng));
                auto child = evaluate(y);
                if (child.obj[party] >= source[c.best].obj[party] - 1e-9 || k < 3)
                    offspring.push_back(child);
            }
        }

        m_party_pop[party].insert(m_party_pop[party].end(), offspring.begin(), offspring.end());
        truncateByFitnessAndDiversity(m_party_pop[party], party, m_config.party_pop_size, 0.018);
    }

    void updateCoordinator(const std::array<std::vector<Cluster>, 2>& clusters,
        const std::array<std::vector<Individual>, 2>& sources) {
        std::vector<Individual> seeds;
        for (int party = 0; party < 2; ++party) {
            for (const auto& c : clusters[party]) {
                if (c.best >= 0 && c.best < static_cast<int>(sources[party].size()))
                    seeds.push_back(sources[party][c.best]);
            }
        }
        for (const auto& a : m_consensus_archive) seeds.push_back(a);
        for (const auto& s : seeds) m_coord_pop.push_back(s);
        truncateByConsensusAndDiversity(m_coord_pop, m_config.coord_pop_size, 0.015);
    }

    void evolveCoordinator() {
        if (m_coord_pop.empty()) return;
        std::uniform_int_distribution<int> I(0, static_cast<int>(m_coord_pop.size()) - 1);
        std::uniform_real_distribution<Real> U(0, 1);
        std::normal_distribution<Real> N(0, 1);
        std::vector<Individual> children;
        children.reserve(m_config.coord_children);

        for (int k = 0; k < m_config.coord_children && remainingEvaluations() > 0; ++k) {
            const auto& base = m_coord_pop[I(m_rng)];
            std::vector<Real> y = base.x;
            if (U(m_rng) < 0.55 && m_coord_pop.size() >= 4) {
                const auto& a = m_coord_pop[I(m_rng)].x;
                const auto& b = m_coord_pop[I(m_rng)].x;
                const auto& c = m_coord_pop[I(m_rng)].x;
                const Real F = 0.55;
                for (size_t d = 0; d < m_dim; ++d)
                    if (U(m_rng) < 0.9) y[d] = clamp01(a[d] + F * (b[d] - c[d]));
            }
            else {
                Real sigma = 0.035 * (1.0 - 0.6 * m_coord_decay);
                for (size_t d = 0; d < m_dim; ++d) y[d] = clamp01(y[d] + sigma * N(m_rng));
            }
            children.push_back(evaluate(y));
        }
        m_coord_decay = std::min<Real>(1.0, m_coord_decay + 1.0 / std::max<int>(1, m_config.iterations));
        m_coord_pop.insert(m_coord_pop.end(), children.begin(), children.end());
        truncateByConsensusAndDiversity(m_coord_pop, m_config.coord_pop_size, 0.012);
    }

    void updateArchives() {
        std::vector<Individual> pool = m_coord_pop;
        for (int p = 0; p < 2; ++p)
            pool.insert(pool.end(), m_party_pop[p].begin(), m_party_pop[p].end());
        for (const auto& ind : pool) addToHistoryArchive(ind);

        std::sort(pool.begin(), pool.end(), [](const Individual& a, const Individual& b) {
            return a.consensus > b.consensus;
        });
        for (const auto& ind : pool) {
            if (ind.consensus < 65.0) continue;
            bool near = false;
            for (auto& a : m_consensus_archive) {
                if (dist(a.x, ind.x) < m_config.archive_radius) {
                    near = true;
                    if (ind.consensus > a.consensus) a = ind;
                    break;
                }
            }
            if (!near && static_cast<int>(m_consensus_archive.size()) < m_config.consensus_archive_max)
                m_consensus_archive.push_back(ind);
        }
        std::sort(m_consensus_archive.begin(), m_consensus_archive.end(), [](const Individual& a, const Individual& b) {
            return a.consensus > b.consensus;
        });
    }

    void addToHistoryArchive(const Individual& ind) {
        const auto key = historyKey(ind);
        if (m_history_keys.find(key) != m_history_keys.end()) return;
        m_history_keys.insert(key);
        m_archive.push_back(ind);
    }

    std::string historyKey(const Individual& ind) const {
        std::ostringstream os;
        os << std::fixed << std::setprecision(12);
        for (Real v : ind.x) os << clamp01(v) << ',';
        return os.str();
    }

    void truncateByFitnessAndDiversity(std::vector<Individual>& pop, int party, int cap, Real min_dist) const {
        std::sort(pop.begin(), pop.end(), [&](const Individual& a, const Individual& b) {
            return a.obj[party] > b.obj[party];
        });
        diverseTruncate(pop, cap, min_dist);
    }

    void truncateByConsensusAndDiversity(std::vector<Individual>& pop, int cap, Real min_dist) const {
        std::sort(pop.begin(), pop.end(), [](const Individual& a, const Individual& b) {
            return a.consensus > b.consensus;
        });
        diverseTruncate(pop, cap, min_dist);
    }

    void diverseTruncate(std::vector<Individual>& pop, int cap, Real min_dist) const {
        std::vector<Individual> kept;
        for (const auto& ind : pop) {
            bool near = false;
            for (const auto& k : kept) {
                if (dist(k.x, ind.x) < min_dist) { near = true; break; }
            }
            if (!near || static_cast<int>(kept.size()) < cap / 2)
                kept.push_back(ind);
            if (static_cast<int>(kept.size()) >= cap) break;
        }
        pop.swap(kept);
    }

    Real bestConsensus() const {
        Real best = 0;
        for (const auto& ind : m_archive) best = std::max(best, ind.consensus);
        for (const auto& ind : m_consensus_archive) best = std::max(best, ind.consensus);
        for (const auto& ind : m_coord_pop) best = std::max(best, ind.consensus);
        return best;
    }

    static Real plotX(Real x) { return x * 7.0 - 3.5; }

    void savePopulation(const std::string& path, const std::vector<Individual>& pop) const {
        std::ofstream out(path);
        out << std::fixed << std::setprecision(10);
        out << "# x0_norm x1_norm x0_plot x1_plot f_party0 f_party1 consensus";
        for (size_t d = 0; d < m_dim; ++d) out << " x" << d << "_norm";
        out << "\n";
        for (const auto& ind : pop) {
            out << ind.x[0] << ' ' << ind.x[1] << ' '
                << plotX(ind.x[0]) << ' ' << plotX(ind.x[1]) << ' '
                << ind.obj[0] << ' ' << ind.obj[1] << ' ' << ind.consensus;
            for (const auto v : ind.x) out << ' ' << v;
            out << "\n";
        }
    }

    void saveNbdEdges(const std::string& path, const std::vector<Individual>& pop, int party) const {
        auto edges = nearestBetterEdges(pop, party);
        std::ofstream out(path);
        out << std::fixed << std::setprecision(10);
        out << "# from_x0 from_x1 to_x0 to_x1 length\n";
        for (const auto& e : edges) {
            out << plotX(pop[e.from].x[0]) << ' ' << plotX(pop[e.from].x[1]) << ' '
                << plotX(pop[e.to].x[0]) << ' ' << plotX(pop[e.to].x[1]) << ' '
                << e.length << "\n";
        }
    }

    void saveRunInfo(const std::string& path) const {
        std::ofstream out(path);
        out << "suite " << m_suite_id << "\n";
        out << "dimension " << m_dim << "\n";
        out << "evaluations " << m_evaluations << "\n";
        out << "max_evaluations " << m_config.max_evaluations << "\n";
        out << "iterations " << m_config.iterations << "\n";
        out << "party_pop_size " << m_config.party_pop_size << "\n";
        out << "coord_pop_size " << m_config.coord_pop_size << "\n";
        out << "history_archive_size " << m_archive.size() << "\n";
        out << "consensus_archive_size " << m_consensus_archive.size() << "\n";
        out << "best_consensus " << std::fixed << std::setprecision(10) << bestConsensus() << "\n";
    }

private:
    FreePeaksMultiParty* m_problem = nullptr;
    Environment* m_env = nullptr;
    int m_suite_id = 1;
    std::string m_output_dir;
    size_t m_dim = 2;
    mutable std::mt19937 m_rng;
    int m_evaluations = 0;

    std::array<std::vector<Individual>, 2> m_party_pop;
    std::vector<Individual> m_coord_pop;
    std::vector<Individual> m_archive;
    std::vector<Individual> m_consensus_archive;
    std::unordered_set<std::string> m_history_keys;

    Config m_config;
    Real m_coord_decay = 0;
};

inline void runMultiPartyNbnCmaesSuite(int suite_id, size_t dimension = 2,
    const std::string& output_root = "Visualization/multiparty_nbn_cmaes",
    const Config& config = Config()) {
    registerInstance();

    std::shared_ptr<Environment> env(generateEnvironmentByFactory("free_peaks_multiparty"));
    env->setProblem(generateProblemByFactory("free_peaks_multiparty"));
    ParameterMap pm;
    pm["suite_id"] = suite_id;
    pm["problem_dimension"] = dimension;
    env->problem()->inputParameters().input(pm);
    env->initializeProblem(0.5);

    auto* problem = dynamic_cast<FreePeaksMultiParty*>(env->problem());
    if (!problem) throw std::runtime_error("free_peaks_multiparty initialization failed");

    const auto& spec = problem->currentSpec();
    std::filesystem::path out_dir = std::filesystem::path(output_root) /
        ("D" + std::to_string(dimension)) /
        ("suite_" + std::to_string(suite_id) + "_" + spec.name);
    Solver solver(problem, env.get(), suite_id, out_dir.string(), config);
    solver.run();
}

inline void runMultiPartyNbnCmaes(int suite_id = 0, size_t dimension = 2,
    const std::string& output_root = "Visualization/multiparty_nbn_cmaes",
    const Config& config = Config()) {
    if (suite_id >= 1 && suite_id <= 12) {
        runMultiPartyNbnCmaesSuite(suite_id, dimension, output_root, config);
        return;
    }
    for (int id = 1; id <= 12; ++id) runMultiPartyNbnCmaesSuite(id, dimension, output_root, config);
}

inline void runMultiPartyNbnCmaesAllDimensions(
    const std::string& output_root = "Visualization/multiparty_nbn_cmaes",
    const Config& config = Config()) {
    for (size_t dim : { size_t(2), size_t(3), size_t(5), size_t(10) }) {
        std::cout << ">>> [MultiParty-NBN-CMAES] Dimension " << dim << "\n";
        runMultiPartyNbnCmaes(0, dim, output_root, config);
    }
}

} // namespace multiparty_nbn_cmaes
} // namespace ofec

#endif
