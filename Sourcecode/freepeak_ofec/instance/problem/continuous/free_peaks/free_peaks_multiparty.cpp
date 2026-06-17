#include "free_peaks_multiparty.h"
#include "factory.h"
#include "subproblem/subproblem.h"
#include "subproblem/distance/distance_base.h"
#include "subproblem/function/one_peak_function.h"
#include "subproblem/function/one_peak/one_peak_base.h"
#include "subproblem/transform/transform_x/transform_x_base.h"
#include "subproblem/transform/transform_y/transform_y_base.h"
#include "../../../../core/environment/environment.h"
#include <algorithm>
#include <cmath>
#include <memory>
#include <limits>
#include <numeric>
#include <sstream>

namespace ofec {

    namespace {
        FreePeaksMultiParty::PartyPeakSpec pk(
            const std::string& name, Real x, Real y, Real h,
            const std::string& shape, Real cond, Real bias, Real rot,
            Real linkage = 0, Real asym = 0, Real radius = 32)
        {
            FreePeaksMultiParty::PartyPeakSpec p;
            p.name = name;
            p.center = { (x + 3.5) / 7.0, (y + 3.5) / 7.0 };
            p.height = h;
            p.shape = shape;
            p.condition = cond;
            p.bias = bias;
            p.rotation = rot;
            p.linkage = linkage;
            p.asymmetry = asym;
            p.radius = radius;
            return p;
        }

        std::vector<std::pair<std::string, std::vector<std::pair<std::string, double>>>> tree(
            std::initializer_list<std::pair<std::string, double>> leaves)
        {
            std::vector<std::pair<std::string, std::vector<std::pair<std::string, double>>>> t;
            t.push_back({ "root", std::vector<std::pair<std::string, double>>(leaves) });
            return t;
        }

        Real intersectionVolume(const std::vector<std::pair<Real, Real>>& a,
            const std::vector<std::pair<Real, Real>>& b)
        {
            Real volume = 1;
            for (size_t d = 0; d < a.size(); ++d) {
                const Real lo = std::max(a[d].first, b[d].first);
                const Real hi = std::min(a[d].second, b[d].second);
                if (hi <= lo) return 0;
                volume *= hi - lo;
            }
            return volume;
        }

        std::vector<Real> sampleInBox(const std::vector<std::pair<Real, Real>>& box,
            Uniform& uniform, Real margin = 0.12)
        {
            std::vector<Real> x(box.size(), 0.5);
            for (size_t d = 0; d < box.size(); ++d) {
                const Real width = box[d].second - box[d].first;
                Real lo = box[d].first + margin * width;
                Real hi = box[d].second - margin * width;
                if (hi <= lo) {
                    lo = box[d].first;
                    hi = box[d].second;
                }
                x[d] = uniform.nextNonStd<Real>(lo, hi);
            }
            return x;
        }

        Real euclideanDistance(const std::vector<Real>& a, const std::vector<Real>& b)
        {
            Real value = 0;
            const size_t n = std::min(a.size(), b.size());
            for (size_t i = 0; i < n; ++i) {
                const Real diff = a[i] - b[i];
                value += diff * diff;
            }
            return std::sqrt(value);
        }

        std::vector<std::pair<Real, Real>> intersectionBox(
            const std::vector<std::pair<Real, Real>>& a,
            const std::vector<std::pair<Real, Real>>& b)
        {
            std::vector<std::pair<Real, Real>> box(a.size());
            for (size_t d = 0; d < a.size(); ++d) {
                box[d] = { std::max(a[d].first, b[d].first), std::min(a[d].second, b[d].second) };
            }
            return box;
        }
    }

    void FreePeaksMultiParty::addInputParameters() {
        Continuous::addInputParameters();
        m_input_parameters.add("suite_id", new RangedInt(m_suite_id, 1, 12, 1));
        m_input_parameters.add("problem_dimension", new RangedSizeT(m_problem_dimension, 2, 100, 2));
    }

    std::array<Real, 2> FreePeaksMultiParty::point(Real x, Real y) { return { (x + 3.5) / 7.0, (y + 3.5) / 7.0 }; }

    const std::vector<FreePeaksMultiParty::SuiteSpec>& FreePeaksMultiParty::suite() {
        static const std::vector<SuiteSpec> specs = {
            { "p1_balanced", "balanced random shared/private optima", {}, {} },
            { "p2_unequal_basins", "random optima with unequal KD-tree basin volumes", {}, {} },
            { "p3_local_conflict", "random shared optima plus party-private local conflicts", {}, {} },
            { "p4_rugged", "random shared optima in rugged/asymmetric transformed subspaces", {}, {} },
            { "p5_rotated", "random optima with rotated and ill-conditioned basins", {}, {} },
            { "p6_deceptive", "random narrow shared optima embedded in broad private basins", {}, {} },
            { "p7_hierarchical", "random broad shared basins plus sharp local peaks", {}, {} },
            { "p8_disconnected", "random common optima across disconnected party partitions", {}, {} },
            { "p9_twenty_shared_optima", "dense random multimodal problem with 20 shared optima", {}, {} },
            { "p10_rotated_rugged_deceptive", "mixed rotated, ill-conditioned, rugged and deceptive shared optima", {}, {} },
            { "p11_hierarchical_disconnected", "mixed hierarchical, disconnected and ill-conditioned shared optima", {}, {} },
            { "p12_full_mixed_multimodal", "dense full-mixed problem with rugged, deceptive, rotated and hierarchical basins", {}, {} }
        };
        return specs;
    }

    std::string FreePeaksMultiParty::partyRegionName(size_t party_id, const std::vector<Real>& x) const {
        return m_party_trees.at(party_id).tree->getRegionName(x);
    }
    const FreePeaksMultiParty::SuiteSpec& FreePeaksMultiParty::currentSpec() const {
        if (!m_current_spec.name.empty()) return m_current_spec;
        const int max_suite = static_cast<int>(suite().size());
        const int id = std::max(1, std::min(max_suite, m_suite_id));
        return suite()[static_cast<size_t>(id - 1)];
    }

    FreePeaksMultiParty::SuiteSpec FreePeaksMultiParty::makeRandomSuiteSpec() {
        const int max_suite = static_cast<int>(suite().size());
        const int id = std::max(1, std::min(max_suite, m_suite_id));
        const auto& meta = suite()[static_cast<size_t>(id - 1)];
        SuiteSpec spec;
        spec.name = meta.name;
        spec.feature = meta.feature;

        const std::array<std::array<size_t, 2>, 8> base_leaf_counts = { {
            { 4, 6 }, { 4, 7 }, { 5, 7 }, { 5, 8 },
            { 4, 7 }, { 5, 8 }, { 5, 9 }, { 6, 10 }
        } };
        const size_t mixed_leaves = id == 9 ? 26 : (id == 10 ? 12 : (id == 11 ? 14 : (id == 12 ? 22 : 7)));
        auto& uniform = random()->uniform;

        for (size_t party = 0; party < 2; ++party) {
            const size_t base_leaves = id <= 8 ?
                base_leaf_counts[static_cast<size_t>(id - 1)][party] : mixed_leaves;
            const size_t leaves_per_party = std::max<size_t>(
                base_leaves, m_problem_dimension + 2 + party * (id <= 8 ? 2 : 0));
            std::vector<std::pair<std::string, double>> leaves;
            leaves.reserve(leaves_per_party);
            for (size_t i = 0; i < leaves_per_party; ++i) {
                std::ostringstream name;
                name << "p" << party << "_r" << i;
                Real weight = 0.6 + uniform.next();
                if (id == 1) weight = party == 0 ? (0.9 + 0.5 * uniform.next()) : (0.35 + 2.20 * uniform.next());
                if (id == 2) weight = (i % 3 == party ? 5.5 + 2.5 * uniform.next() : 0.25 + 0.65 * uniform.next());
                if (id == 3) weight = (i % 2 == party ? 3.2 + 1.7 * uniform.next() : 0.45 + 0.85 * uniform.next());
                if (id == 4) weight = (i % 4 <= 1 ? 0.35 + 0.45 * uniform.next() : 2.8 + 1.8 * uniform.next());
                if (id == 5) weight = (party == 0 ? 0.55 + 1.8 * uniform.next() : 0.25 + 3.4 * uniform.next());
                if (id == 6) weight = (i < 2 + party ? 5.0 + 2.5 * uniform.next() : 0.30 + 0.80 * uniform.next());
                if (id == 7) weight = (i % 4 == 0 ? 7.0 + 2.0 * uniform.next() : 0.30 + 0.90 * uniform.next());
                if (id == 8) weight = (i % 3 == party ? 4.8 + 2.4 * uniform.next() : 0.25 + 0.75 * uniform.next());
                if (id == 10 && i % 4 == party) weight = 3.1;
                if (id == 11 && i < 3) weight = 4.0 - 0.3 * i;
                if (id == 12 && i % 5 < 2) weight = 2.8 + 0.5 * uniform.next();
                leaves.emplace_back(name.str(), static_cast<double>(weight));
            }
            uniform.shuffle(leaves.begin(), leaves.end());
            spec.parties[party].tree.push_back({ "root", leaves });
        }
        return spec;
    }

    void FreePeaksMultiParty::generateRandomPeaks(SuiteSpec& spec) {
        struct Candidate {
            size_t p0 = 0;
            size_t p1 = 0;
            Real volume = 0;
        };

        const int max_suite = static_cast<int>(suite().size());
        const int id = std::max(1, std::min(max_suite, m_suite_id));
        const std::array<size_t, 12> shared_counts = { 2, 2, 1, 2, 1, 2, 2, 2, 20, 6, 8, 16 };
        const size_t desired_shared = shared_counts[static_cast<size_t>(id - 1)];
        auto& uniform = random()->uniform;
        auto unit = [&uniform]() { return uniform.next(); };
        auto tuneForVisibleBasins = [](PartyPeakSpec& peak, int suite_id) {
            Real min_radius = 86;
            Real max_condition = 14;
            Real max_rotation = 0.95;
            Real max_linkage = 14;
            Real max_asymmetry = 1.35;
            Real max_deceptive = 0.32;

            if (suite_id == 1) { min_radius = 76; max_condition = 10; max_rotation = 0.70; }
            if (suite_id == 2) { min_radius = 110; max_condition = 10; max_rotation = 0.65; }
            if (suite_id == 3) { min_radius = 160; max_condition = 5; max_rotation = 0.45; max_deceptive = 0.14; }
            if (suite_id == 4) { min_radius = 92; max_condition = 12; max_linkage = 10; max_asymmetry = 1.10; max_deceptive = 0.24; }
            if (suite_id == 5) { min_radius = 185; max_condition = 5; max_rotation = 0.50; max_deceptive = 0.12; }
            if (suite_id == 6) { min_radius = 75; max_condition = 5; max_rotation = 0.45; max_deceptive = 0.32; }
            if (suite_id == 7) { min_radius = 190; max_condition = 4.5; max_rotation = 0.42; max_deceptive = 0.12; }
            if (suite_id == 8) { min_radius = 180; max_condition = 5.5; max_rotation = 0.48; max_deceptive = 0.14; }
            if (suite_id == 9) { min_radius = 125; max_condition = 7; max_rotation = 0.48; max_linkage = 7; max_asymmetry = 0.85; max_deceptive = 0.16; }
            if (suite_id >= 10) { min_radius = 185; max_condition = 5.5; max_rotation = 0.42; max_linkage = 6; max_asymmetry = 0.80; max_deceptive = 0.14; }
            if (suite_id == 10) { min_radius = 210; max_condition = 4.5; max_rotation = 0.34; max_linkage = 5; max_asymmetry = 0.70; max_deceptive = 0.12; }
            if (suite_id == 11) { min_radius = 195; max_condition = 5.0; max_rotation = 0.38; max_linkage = 5.5; max_asymmetry = 0.75; max_deceptive = 0.12; }
            if (suite_id == 12) { min_radius = 220; max_condition = 4.5; max_rotation = 0.34; max_linkage = 5; max_asymmetry = 0.70; max_deceptive = 0.12; }

            peak.radius = std::max(peak.radius, min_radius);
            peak.condition = std::max<Real>(1.0, std::min(peak.condition, max_condition));
            if (peak.rotation != 0) {
                const Real sign = peak.rotation < 0 ? -1.0 : 1.0;
                peak.rotation = sign * std::min(std::abs(peak.rotation), max_rotation);
            }
            peak.linkage = std::min(peak.linkage, max_linkage);
            peak.asymmetry = std::min(peak.asymmetry, max_asymmetry);
            peak.deceptive = std::min(peak.deceptive, max_deceptive);
        };

        std::array<std::vector<std::string>, 2> leaf_names;
        for (size_t party = 0; party < 2; ++party) {
            for (const auto& leaf : spec.parties[party].tree.front().second)
                leaf_names[party].push_back(leaf.first);
        }

        std::vector<Candidate> candidates;
        for (size_t i = 0; i < leaf_names[0].size(); ++i) {
            const auto& box0 = m_party_trees[0].tree->getBox(leaf_names[0][i]);
            for (size_t j = 0; j < leaf_names[1].size(); ++j) {
                const auto& box1 = m_party_trees[1].tree->getBox(leaf_names[1][j]);
                const Real vol = intersectionVolume(box0, box1);
                if (vol > 1e-12) candidates.push_back({ i, j, vol });
            }
        }
        std::sort(candidates.begin(), candidates.end(),
            [](const Candidate& a, const Candidate& b) { return a.volume > b.volume; });

        std::vector<int> used0(leaf_names[0].size(), 0), used1(leaf_names[1].size(), 0);
        size_t shared_count = 0;
        for (const auto& c : candidates) {
            if (shared_count >= desired_shared) break;
            if (used0[c.p0] || used1[c.p1]) continue;
            used0[c.p0] = 1;
            used1[c.p1] = 1;

            const auto& box0 = m_party_trees[0].tree->getBox(leaf_names[0][c.p0]);
            const auto& box1 = m_party_trees[1].tree->getBox(leaf_names[1][c.p1]);
            auto shared_box = intersectionBox(box0, box1);
            auto x = sampleInBox(shared_box, uniform, 0.18);
            spec.shared_optima.push_back(x);

            for (size_t party = 0; party < 2; ++party) {
                PartyPeakSpec peak;
                peak.name = party == 0 ? leaf_names[0][c.p0] : leaf_names[1][c.p1];
                peak.center = x;
                peak.shape = "mp";
                peak.height = 90;
                peak.radius = 48 + 36 * unit();
                if (id == 1) {
                    const size_t mode = (shared_count + party) % 3;
                    peak.radius = (party == 0 ? 34 + 26 * mode : 96 - 20 * mode) + 22 * unit();
                    peak.condition = 1.0 + (party == 0 ? mode : 2 - mode) * 1.6 + 2.5 * unit();
                    peak.height = 88 + 3 * unit();
                }
                if (id == 2) {
                    const size_t mode = shared_count % 4;
                    if (party == 0) {
                        peak.radius = (mode == 0 ? 18 : (mode == 1 ? 180 : 44)) + 20 * unit();
                        peak.condition = 2 + 12 * unit();
                    }
                    else {
                        peak.radius = (mode == 0 ? 170 : (mode == 1 ? 22 : 120)) + 30 * unit();
                        peak.condition = 8 + 24 * unit();
                    }
                }
                if (id == 3) {
                    const size_t mode = (shared_count * 2 + party) % 5;
                    peak.radius = (mode < 2 ? 64 : 112) + (party == 0 ? 34 : 42) * unit();
                    peak.condition = (party == 0 ? 4 : 8) + mode * 2.5 + 10 * unit();
                    peak.rotation = (party == 0 ? 1 : -1) * (0.25 + 0.80 * unit());
                    if (party == 1 || mode >= 3) peak.deceptive = 0.06 + 0.14 * unit();
                }
                if (id == 4) {
                    const size_t mode = (shared_count + 2 * party) % 4;
                    peak.radius = (party == 0 ? 30 + 18 * mode : 92 - 14 * mode) + 26 * unit();
                    peak.linkage = (party == 0 ? 10.0 : 18.0) + mode * 6.0 + 20.0 * unit();
                    peak.asymmetry = (party == 0 ? 1.00 : 1.45) + 0.55 * mode + 1.20 * unit();
                    peak.deceptive = (mode % 2 == 0 ? 0.18 : 0.34) + (party == 0 ? 0.22 : 0.34) * unit();
                }
                if (id == 5 || id == 8) {
                    const size_t mode = (shared_count * 3 + party) % 5;
                    if (id == 5) {
                        peak.radius = (party == 0 ? 72 + 12 * mode : 118 - 8 * mode) + 36 * unit();
                        peak.condition = (party == 0 ? 10 : 16) + mode * 4 + 16 * unit();
                        peak.rotation = (party == 0 ? 1 : -1) * (0.55 + 0.28 * mode + 1.05 * unit());
                    }
                    else {
                        peak.radius = (mode == 0 || mode == 3 ? 76 : 150) + (party == 0 ? 38 : 52) * unit();
                        peak.condition = (party == 0 ? 8 : 14) + mode * 3 + 18 * unit();
                        peak.rotation = (party == 0 ? 1 : -1) * (0.30 + 0.20 * mode + 0.85 * unit());
                        peak.deceptive = (mode % 2 == 0 ? 0.08 : 0.18) + (party == 0 ? 0.14 : 0.22) * unit();
                    }
                }
                if (id == 6) {
                    const size_t mode = (shared_count + party) % 4;
                    peak.radius = (mode == 0 ? 14 : (mode == 1 ? 150 : 30)) + (party == 0 ? 14 : 30) * unit();
                    peak.deceptive = (party == 0 ? 0.40 : 0.55) + mode * 0.08 + 0.24 * unit();
                    peak.condition = (party == 0 ? 14 : 30) + mode * 8 + 42 * unit();
                }
                if (id == 7) {
                    const size_t mode = shared_count % 4;
                    if (party == 0) {
                        peak.radius = (mode == 0 ? 175 : (mode == 1 ? 58 : 118)) + 34 * unit();
                        peak.condition = (mode == 0 ? 2 : 9) + 16 * unit();
                    }
                    else {
                        peak.radius = (mode == 0 ? 62 : (mode == 1 ? 170 : 92)) + 42 * unit();
                        peak.condition = (mode == 1 ? 3 : 12) + 18 * unit();
                    }
                }
                if (id == 9) {
                    const size_t mode = shared_count % 5;
                    if (mode == 0) {
                        peak.radius = 34 + 24 * unit();
                        peak.condition = 22 + 58 * unit();
                        peak.rotation = (party == 0 ? 1 : -1) * (0.85 + 1.55 * unit());
                    }
                    else if (mode == 1) {
                        peak.radius = 45 + 35 * unit();
                        peak.linkage = 6.0 + 13.0 * unit();
                        peak.asymmetry = 0.65 + 0.9 * unit();
                        peak.deceptive = 0.12 + 0.18 * unit();
                    }
                    else if (mode == 2) {
                        peak.deceptive = 0.24 + 0.28 * unit();
                        peak.radius = 24 + 16 * unit();
                        peak.condition = 7 + 18 * unit();
                    }
                    else if (mode == 3) {
                        peak.radius = 95 + 45 * unit();
                        peak.condition = 4 + 10 * unit();
                    }
                    else {
                        peak.radius = 35 + 28 * unit();
                        peak.condition = 18 + 42 * unit();
                        peak.rotation = (party == 0 ? -1 : 1) * (0.65 + 1.35 * unit());
                        peak.deceptive = 0.16 + 0.24 * unit();
                    }
                }
                if (id == 10) {
                    const size_t mode = shared_count % 6;
                    if (party == 0) {
                        peak.radius = (mode % 2 == 0 ? 115 : 55) + 35 * unit();
                        peak.condition = (mode < 3 ? 8 : 26) + 22 * unit();
                        peak.rotation = 0.85 + 1.45 * unit();
                        if (mode == 0 || mode == 3) {
                            peak.linkage = 8.0 + 12.0 * unit();
                            peak.asymmetry = 0.60 + 0.95 * unit();
                        }
                        if (mode == 1 || mode == 4) peak.deceptive = 0.24 + 0.32 * unit();
                        if (mode == 5) {
                            peak.linkage = 3.0 + 7.0 * unit();
                            peak.deceptive = 0.14 + 0.24 * unit();
                        }
                    }
                    else {
                        peak.radius = (mode % 2 == 0 ? 38 : 145) + 28 * unit();
                        peak.condition = (mode < 3 ? 34 : 10) + 30 * unit();
                        peak.rotation = -(0.25 + 1.95 * unit());
                        if (mode == 0 || mode == 2) peak.deceptive = 0.30 + 0.30 * unit();
                        if (mode == 1 || mode == 4 || mode == 5) {
                            peak.linkage = 10.0 + 16.0 * unit();
                            peak.asymmetry = 0.85 + 1.20 * unit();
                        }
                    }
                }
                if (id == 11) {
                    const size_t mode = shared_count % 8;
                    if (party == 0) {
                        peak.radius = (mode % 4 == 0 ? 150 : (mode % 4 == 1 ? 28 : 72)) + 22 * unit();
                        peak.condition = (mode % 4 == 0 ? 3 : 18) + 34 * unit();
                        peak.rotation = -(0.20 + 1.25 * unit());
                        if (mode % 4 == 1 || mode % 4 == 3) peak.deceptive = 0.16 + 0.30 * unit();
                        if (mode % 4 >= 2) {
                            peak.linkage = 4.0 + 11.0 * unit();
                            peak.asymmetry = 0.25 + 0.95 * unit();
                        }
                    }
                    else {
                        peak.radius = (mode % 4 == 0 ? 30 : (mode % 4 == 1 ? 135 : 42)) + 24 * unit();
                        peak.condition = (mode % 4 == 1 ? 4 : 30) + 42 * unit();
                        peak.rotation = 0.55 + 1.60 * unit();
                        if (mode % 2 == 0) peak.deceptive = 0.20 + 0.34 * unit();
                        if (mode % 4 == 0 || mode % 4 == 3) {
                            peak.linkage = 9.0 + 17.0 * unit();
                            peak.asymmetry = 0.70 + 1.20 * unit();
                        }
                    }
                }
                if (id == 12) {
                    const size_t mode = shared_count % 8;
                    if (party == 0) {
                        peak.radius = (mode == 0 || mode == 5 ? 140 : (mode == 3 ? 30 : 65)) + 35 * unit();
                        peak.condition = (mode == 1 || mode == 6 ? 40 : 8) + 32 * unit();
                        peak.rotation = 0.15 + 2.20 * unit();
                        if (mode == 1 || mode == 3 || mode == 7) peak.deceptive = 0.24 + 0.36 * unit();
                        if (mode == 2 || mode == 4 || mode == 7) {
                            peak.linkage = 5.0 + 18.0 * unit();
                            peak.asymmetry = 0.45 + 1.35 * unit();
                        }
                    }
                    else {
                        peak.radius = (mode == 0 || mode == 5 ? 36 : (mode == 3 ? 120 : 48)) + 38 * unit();
                        peak.condition = (mode == 2 || mode == 7 ? 46 : 12) + 38 * unit();
                        peak.rotation = -(0.55 + 1.90 * unit());
                        if (mode == 0 || mode == 4 || mode == 6) peak.deceptive = 0.20 + 0.42 * unit();
                        if (mode == 1 || mode == 5 || mode == 6) {
                            peak.linkage = 8.0 + 20.0 * unit();
                            peak.asymmetry = 0.70 + 1.40 * unit();
                        }
                    }
                }
                tuneForVisibleBasins(peak, id);
                spec.parties[party].peaks.push_back(peak);
            }
            ++shared_count;
        }

        std::array<int, 2> p6_neighbor_deceptive_done { 0, 0 };
        for (size_t party = 0; party < 2; ++party) {
            for (size_t i = 0; i < leaf_names[party].size(); ++i) {
                if ((party == 0 && used0[i]) || (party == 1 && used1[i])) continue;
                const auto& box = m_party_trees[party].tree->getBox(leaf_names[party][i]);
                PartyPeakSpec peak;
                peak.name = leaf_names[party][i];
                peak.center = sampleInBox(box, uniform, 0.14);
                peak.shape = "mp";
                peak.height = 62 + 24 * unit();
                peak.radius = 40 + 58 * unit();
                if (id == 1) {
                    const size_t mode = (i + 2 * party) % 4;
                    peak.height = 64 + (party == 0 ? 18 : 24) * unit();
                    peak.radius = (party == 0 ? 28 + 24 * mode : 108 - 18 * mode) + 26 * unit();
                    peak.condition = 1.0 + mode * (party == 0 ? 1.2 : 1.8) + 2.5 * unit();
                }
                if (id == 2) {
                    const size_t mode = (i + party) % 5;
                    peak.height = party == 0 ? 58 + 26 * unit() : 64 + 20 * unit();
                    if (party == 0) {
                        peak.radius = (mode == 0 ? 185 : (mode == 2 ? 16 : 42)) + 32 * unit();
                        peak.condition = 2 + 16 * unit();
                    }
                    else {
                        peak.radius = (mode == 0 ? 24 : (mode == 2 ? 190 : 115)) + 44 * unit();
                        peak.condition = 8 + 30 * unit();
                    }
                }
                if (id == 3) {
                    const size_t mode = (i * 2 + party) % 5;
                    peak.height = party == 0 ? 66 + 18 * unit() : 56 + 30 * unit();
                    peak.radius = (mode < 2 ? 18 : 105) + (party == 0 ? 22 : 38) * unit();
                    peak.condition = (party == 0 ? 12 : 26) + mode * 9 + 36 * unit();
                    peak.rotation = (party == 0 ? 1 : -1) * (0.25 + 1.80 * unit());
                    if (mode == 1 || mode == 3 || party == 1) peak.deceptive = 0.10 + 0.28 * unit();
                }
                if (id == 6 && i < 2) {
                    peak.radius = 105 + 30 * unit();
                    peak.height = 68 + 10 * unit();
                }
                if (id == 5 || id == 8) {
                    const size_t mode = (i * 3 + party) % 6;
                    if (id == 5) {
                        peak.radius = (party == 0 ? 34 + 18 * mode : 145 - 16 * mode) + 38 * unit();
                        peak.condition = (party == 0 ? 34 : 72) + mode * 16 + 62 * unit();
                        peak.rotation = (party == 0 ? 1 : -1) * (0.75 + 0.65 * mode + 2.45 * unit());
                    }
                    else {
                        peak.radius = (mode == 0 || mode == 4 ? 22 : 142) + (party == 0 ? 34 : 54) * unit();
                        peak.condition = (party == 0 ? 28 : 58) + mode * 10 + 62 * unit();
                        peak.rotation = (party == 0 ? 1 : -1) * (0.55 + 0.50 * mode + 2.10 * unit());
                        peak.deceptive = (mode % 2 == 0 ? 0.16 : 0.34) + (party == 0 ? 0.28 : 0.40) * unit();
                    }
                }
                if (id == 4) {
                    const size_t mode = (i + 3 * party) % 5;
                    peak.radius = (party == 0 ? 24 + 18 * mode : 106 - 14 * mode) + 30 * unit();
                    peak.linkage = (party == 0 ? 9.0 : 18.0) + mode * 5.5 + 24.0 * unit();
                    peak.asymmetry = (party == 0 ? 0.90 : 1.35) + 0.45 * mode + 1.30 * unit();
                    peak.deceptive = (mode % 2 == 0 ? 0.16 : 0.34) + (party == 0 ? 0.24 : 0.36) * unit();
                }
                if (id == 6) {
                    const size_t mode = (i + 2 * party) % 5;
                    peak.height = party == 0 ? 60 + 22 * unit() : 54 + 30 * unit();
                    peak.radius = (mode == 0 ? 180 : (mode == 2 ? 14 : 38)) + (party == 0 ? 24 : 42) * unit();
                    peak.condition = (party == 0 ? 16 : 34) + mode * 8 + 48 * unit();
                    peak.deceptive = (party == 0 ? 0.42 : 0.58) + mode * 0.07 + 0.28 * unit();
                    if (!p6_neighbor_deceptive_done[party] && !spec.shared_optima.empty()) {
                        const auto& opt = spec.shared_optima.front();
                        std::vector<Real> near_opt(opt.size());
                        for (size_t d = 0; d < opt.size(); ++d) {
                            const Real offset = (d % 2 == party ? 0.10 : -0.10);
                            const Real width = box[d].second - box[d].first;
                            const Real lo = box[d].first + 0.08 * width;
                            const Real hi = box[d].second - 0.08 * width;
                            near_opt[d] = std::max(lo, std::min(hi, opt[d] + offset));
                        }
                        if (euclideanDistance(near_opt, opt) < 0.24) {
                            peak.center = near_opt;
                            peak.height = 80 + 4 * unit();
                            peak.radius = 285 + 45 * unit();
                            peak.condition = 2.0 + 2.0 * unit();
                            peak.rotation = (party == 0 ? 1 : -1) * (0.08 + 0.18 * unit());
                            peak.deceptive = 0.22 + 0.10 * unit();
                            p6_neighbor_deceptive_done[party] = 1;
                        }
                    }
                }
                if (id == 7) {
                    const size_t mode = (i + party) % 4;
                    peak.height = 60 + (party == 0 ? 24 : 18) * unit();
                    if (party == 0) {
                        peak.radius = (mode == 0 ? 205 : (mode == 1 ? 15 : 96)) + 32 * unit();
                        peak.condition = (mode == 0 ? 2 : 42) + 60 * unit();
                    }
                    else {
                        peak.radius = (mode == 0 ? 18 : (mode == 1 ? 195 : 36)) + 42 * unit();
                        peak.condition = (mode == 1 ? 3 : 55) + 65 * unit();
                    }
                }
                if (id == 10) {
                    const size_t mode = (i * 3 + party) % 7;
                    peak.height = party == 0 ? 60 + 24 * unit() : 54 + 32 * unit();
                    if (party == 0) {
                        peak.radius = (mode == 0 || mode == 4 ? 128 : (mode == 2 ? 26 : 58)) + 28 * unit();
                        peak.condition = (mode == 1 || mode == 5 ? 48 : 12) + 38 * unit();
                        peak.rotation = 0.55 + 2.35 * unit();
                        if (mode == 2 || mode == 6) peak.deceptive = 0.28 + 0.34 * unit();
                        if (mode == 0 || mode == 3 || mode == 5) {
                            peak.linkage = 7.0 + 18.0 * unit();
                            peak.asymmetry = 0.70 + 1.25 * unit();
                        }
                    }
                    else {
                        peak.radius = (mode == 0 || mode == 4 ? 30 : (mode == 2 ? 145 : 46)) + 34 * unit();
                        peak.condition = (mode == 3 || mode == 6 ? 58 : 18) + 46 * unit();
                        peak.rotation = -(0.20 + 2.65 * unit());
                        if (mode == 0 || mode == 5) peak.deceptive = 0.34 + 0.34 * unit();
                        if (mode == 1 || mode == 2 || mode == 4) {
                            peak.linkage = 10.0 + 22.0 * unit();
                            peak.asymmetry = 0.85 + 1.45 * unit();
                        }
                    }
                }
                if (id == 11) {
                    const size_t mode = (i + 2 * party) % 8;
                    peak.height = party == 0 ? 57 + 27 * unit() : 62 + 22 * unit();
                    if (party == 0) {
                        peak.radius = (mode < 2 ? 165 : (mode == 3 || mode == 6 ? 24 : 82)) + 30 * unit();
                        peak.condition = (mode < 2 ? 3 : 25) + 55 * unit();
                        peak.rotation = -(0.15 + 1.75 * unit());
                        if (mode == 2 || mode == 5 || mode == 7) peak.deceptive = 0.22 + 0.36 * unit();
                        if (mode == 1 || mode == 4 || mode == 6) {
                            peak.linkage = 5.0 + 20.0 * unit();
                            peak.asymmetry = 0.40 + 1.50 * unit();
                        }
                    }
                    else {
                        peak.radius = (mode < 2 ? 28 : (mode == 3 || mode == 6 ? 155 : 38)) + 36 * unit();
                        peak.condition = (mode == 0 || mode == 4 ? 62 : 9) + 50 * unit();
                        peak.rotation = 0.65 + 2.10 * unit();
                        if (mode == 0 || mode == 3 || mode == 6) peak.deceptive = 0.26 + 0.40 * unit();
                        if (mode == 2 || mode == 5 || mode == 7) {
                            peak.linkage = 11.0 + 24.0 * unit();
                            peak.asymmetry = 0.75 + 1.55 * unit();
                        }
                    }
                }
                if (id == 12) {
                    const size_t mode = (i * 5 + party * 3) % 10;
                    peak.height = party == 0 ? 54 + 34 * unit() : 50 + 38 * unit();
                    if (party == 0) {
                        peak.radius = (mode == 0 || mode == 7 ? 172 : (mode == 2 || mode == 5 ? 18 : 55)) + 32 * unit();
                        peak.condition = (mode == 1 || mode == 8 ? 70 : 10) + 58 * unit();
                        peak.rotation = (mode % 2 == 0 ? 1 : -1) * (0.45 + 2.85 * unit());
                        if (mode == 3 || mode == 5 || mode == 9) peak.deceptive = 0.32 + 0.42 * unit();
                        if (mode == 0 || mode == 4 || mode == 6 || mode == 8) {
                            peak.linkage = 8.0 + 28.0 * unit();
                            peak.asymmetry = 0.80 + 1.70 * unit();
                        }
                    }
                    else {
                        peak.radius = (mode == 0 || mode == 7 ? 22 : (mode == 2 || mode == 5 ? 160 : 44)) + 40 * unit();
                        peak.condition = (mode == 4 || mode == 9 ? 78 : 14) + 62 * unit();
                        peak.rotation = (mode % 2 == 0 ? -1 : 1) * (0.30 + 3.05 * unit());
                        if (mode == 1 || mode == 6 || mode == 8) peak.deceptive = 0.36 + 0.42 * unit();
                        if (mode == 2 || mode == 3 || mode == 5 || mode == 7) {
                            peak.linkage = 12.0 + 30.0 * unit();
                            peak.asymmetry = 0.95 + 1.85 * unit();
                        }
                    }
                }
                tuneForVisibleBasins(peak, id);
                spec.parties[party].peaks.push_back(peak);
            }
        }
    }

    void FreePeaksMultiParty::initialize_(Environment* env) {
        Continuous::initialize_(env);
        FreePeaks::registerFP();
        buildAssignedSuite(env);
    }

    void FreePeaksMultiParty::buildAssignedSuite(Environment* env) {
        m_problem_dimension = std::max<size_t>(2, m_problem_dimension);
        setSizes(m_problem_dimension, 2, 0);
        m_generation_type = "assigned";
        m_file_name = "multiparty/free_peaks_multiparty_suite_" + std::to_string(m_suite_id)
            + "_D" + std::to_string(m_problem_dimension) + ".txt";
        m_party_trees.clear();
		m_party_trees.resize(2);
        m_current_spec = makeRandomSuiteSpec();
        for (size_t p = 0; p < 2; ++p) {
            setPartyKDTree(p, m_current_spec.parties[p].tree);
        }
        generateRandomPeaks(m_current_spec);
        for (size_t p = 0; p < 2; ++p) {
            for (const auto& peak_spec : m_current_spec.parties[p].peaks) setPartySubproblem(p, peak_spec);
        }
        for (auto& party_tree : m_party_trees) {
            for (auto& it : party_tree.name_box_subproblem) {
                if (it.second.second) it.second.second->bindData();
            }
        }
        calculateOptima();
    }

    void FreePeaksMultiParty::setPartyKDTree(size_t party_id,
        const std::vector<std::pair<std::string, std::vector<std::pair<std::string, double>>>>& tree_data)
    {
        std::vector<std::pair<Real, Real>> ranges(numberVariables());
        for (size_t i = 0; i < numberVariables(); ++i) ranges[i] = m_domain[i].limit;
        m_party_trees[party_id].clear();
        m_party_trees[party_id].createTree(ranges, tree_data);
    }

    void FreePeaksMultiParty::setPartySubproblem(size_t party_id, const PartyPeakSpec& peak_spec) {
        using namespace free_peaks;
        ParameterMap spm;
        spm["subspace"] = peak_spec.name;
        spm["generation_type"] = std::string("assigned");
        spm["dataFile1"] = std::string("multiparty/") + peak_spec.name + ".txt";
        auto subpro(Subproblem::create());
        subpro->initialize(spm, this);

        {
            auto dis(FactoryFP<DistanceBase>::produce("flat_border"));
            ParameterMap dp;
            dis->initialize(this, peak_spec.name, dp);
            subpro->setDistance(dis);
        }

        {
            ParameterMap fpm;
            fpm["generation_type"] = std::string("assigned");
            fpm["dataFile1"] = std::string("multiparty/") + peak_spec.name + "_fn.txt";
            auto func(FactoryFP<FunctionBase>::produce("one_peak"));
            func->initialize(this, peak_spec.name, fpm);
            auto* opf = dynamic_cast<OnePeakFunction*>(func);
            for (size_t obj = 0; obj < numberObjectives(); ++obj) {
                auto op(FactoryFP<OnePeakBase>::produce(peak_spec.shape));
                ParameterMap opp;
                opp["center_type"] = std::string("assigned");
                opp["height"] = peak_spec.height;
                opp["radius"] = peak_spec.radius;
                const auto& box = m_party_trees[party_id].tree->getBox(peak_spec.name);
                std::vector<Real> center(numberVariables(), Real(0));
                for (size_t d = 0; d < numberVariables(); ++d) {
                    const Real normalized = d < peak_spec.center.size() ? peak_spec.center[d] : Real(0.5);
                    center[d] = (normalized - box[d].first) / (box[d].second - box[d].first) * 200.0 - 100.0;
                }
                opp["center_postion"] = center;
                op->initialize(this, peak_spec.name, opp);
                opf->addOnePeaks(op);
            }
            subpro->setFunction(func);
        }

        if (peak_spec.bias != 0 || peak_spec.rotation != 0 || peak_spec.condition != 1) {
            auto t(FactoryFP<X_TransformBase>::produce("MapXPartyBias"));
            ParameterMap tp;
            tp["party_id"] = static_cast<int>(party_id);
            tp["magnitude"] = peak_spec.bias;
            tp["rotation_angle"] = peak_spec.rotation;
            tp["condition_number"] = 1.0;
            t->initialize(this, peak_spec.name, tp);
            subpro->addVariableTransform(t);
        }
        if (peak_spec.condition != 1) {
            auto t(FactoryFP<X_TransformBase>::produce("MapXIllConditioning"));
            ParameterMap tp;
            tp["condition"] = peak_spec.condition;
            t->initialize(this, peak_spec.name, tp);
            subpro->addVariableTransform(t);
        }
        if (peak_spec.linkage != 0) {
            auto t(FactoryFP<X_TransformBase>::produce("MapXLinkage"));
            ParameterMap tp;
            tp["beta"] = peak_spec.linkage;
            t->initialize(this, peak_spec.name, tp);
            subpro->addVariableTransform(t);
        }
        if (peak_spec.asymmetry != 0) {
            auto t(FactoryFP<X_TransformBase>::produce("MapXAsymmetry"));
            ParameterMap tp;
            tp["beta"] = peak_spec.asymmetry;
            t->initialize(this, peak_spec.name, tp);
            subpro->addVariableTransform(t);
        }
        if (peak_spec.deceptive != 0) {
            auto t(FactoryFP<X_TransformBase>::produce("MapXDeceptive"));
            ParameterMap tp;
            tp["alpha"] = peak_spec.deceptive;
            t->initialize(this, peak_spec.name, tp);
            subpro->addVariableTransform(t);
        }


        m_party_trees[party_id].name_box_subproblem.at(peak_spec.name).second.reset(subpro);
    }

    free_peaks::Subproblem* FreePeaksMultiParty::subproblem(const std::string& subspace_name) const {
        for (const auto& party_tree : m_party_trees) {
            auto it = party_tree.name_box_subproblem.find(subspace_name);
            if (it != party_tree.name_box_subproblem.end()) return it->second.second.get();
        }
        return FreePeaks::subproblem(subspace_name);
    }

    const std::vector<std::pair<Real, Real>>& FreePeaksMultiParty::subspaceBox(const std::string& subspace_name) const {
        for (const auto& party_tree : m_party_trees) {
            auto it = party_tree.name_box_subproblem.find(subspace_name);
            if (it != party_tree.name_box_subproblem.end()) return party_tree.tree->getBox(subspace_name);
        }
        return FreePeaks::subspaceBox(subspace_name);
    }
    void FreePeaksMultiParty::evaluate(const VariableBase& vars, std::vector<Real>& objs, std::vector<Real>& cons) const {
        auto& x = dynamic_cast<const VariableVector<>&>(vars);
        std::vector<Real> x_vec(x.begin(), x.end());
        objs.assign(numberObjectives(), -std::numeric_limits<Real>::max());
        if (boundaryViolated(vars)) {
            objs = m_default_objective;
            cons = m_default_contrait;
            return;
        }
        for (size_t party = 0; party < numberObjectives(); ++party) {
            auto name = m_party_trees[party].tree->getRegionName(x_vec);
            std::vector<Real> party_objs(numberObjectives(), -std::numeric_limits<Real>::max());
            std::vector<Real> party_cons;
            m_party_trees[party].name_box_subproblem.at(name).second->evaluate(x_vec, party_objs, party_cons);
            objs[party] = party_objs[party];
        }
        cons.clear();
    }

    void FreePeaksMultiParty::calculateOptima() {
        auto new_optima = new Optima<>();
        const auto& spec = currentSpec();
        for (const auto& pt : spec.shared_optima) {
            Solution<> sol(numberObjectives(), numberConstraints(), numberVariables());
            sol.variable().vector() = pt;
            evaluate(sol.variable(), sol.objective(), sol.constraint());
            new_optima->appendSolution(sol);
        }
        m_optima.reset(new_optima);
        m_objective_accuracy = 1e-3;
    }

    void FreePeaksMultiParty::updateOptima(Environment*) {
        calculateOptima();
    }
}
