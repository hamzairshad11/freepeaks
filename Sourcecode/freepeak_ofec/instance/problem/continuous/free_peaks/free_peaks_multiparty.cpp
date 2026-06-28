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
        // Returns m_suite_id clamped to [1, suite().size()].
        // Used by currentSpec(), makeRandomSuiteSpec(), and generateRandomPeaks()
        // to avoid repeating the same two-line clamp expression in each function.
        int clampSuiteId(int suite_id) {
            return std::max(1, std::min(static_cast<int>(FreePeaksMultiParty::suite().size()), suite_id));
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
            { "p3_neighbor_deception", "broad private peaks placed near shared optima creating spatial deception", {}, {} },
            { "p4_rugged", "random shared optima in rugged/asymmetric transformed subspaces", {}, {} },
            { "p5_rotated", "random optima with rotated and ill-conditioned basins", {}, {} },
            { "p6_strong_deception", "narrow shared optima (h=90) embedded in very broad high-fitness private basins (h~87)", {}, {} },
            { "p7_multiscale_hierarchy", "extreme k-d weight contrast: broad private peaks in large leaves, shared peaks in tiny leaves", {}, {} },
            { "p8_party_asymmetric", "party-0 broad smooth basins vs party-1 narrow ill-conditioned basins sharing common optima", {}, {} },
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
        return suite()[static_cast<size_t>(clampSuiteId(m_suite_id) - 1)];
    }

    FreePeaksMultiParty::SuiteSpec FreePeaksMultiParty::makeRandomSuiteSpec() {
        const int id = clampSuiteId(m_suite_id);
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
                // Suite 3 Neighbor Deception: balanced layout, uniform leaf sizes
                if (id == 3) weight = 0.7 + 0.6 * uniform.next();
                if (id == 4) weight = (i % 4 <= 1 ? 0.35 + 0.45 * uniform.next() : 2.8 + 1.8 * uniform.next());
                if (id == 5) weight = (party == 0 ? 0.55 + 1.8 * uniform.next() : 0.25 + 3.4 * uniform.next());
                // Suite 6 Strong Deception: shared-opt leaves (first 2) are TINY, rest are LARGE
                if (id == 6) weight = (i < 2 ? 0.12 + 0.06 * uniform.next() : 4.2 + 2.2 * uniform.next());
                // Suite 7 Multi-scale Hierarchy: extreme contrast (1 in 3 leaves is very large)
                if (id == 7) weight = (i % 3 == 0 ? 9.0 + 3.0 * uniform.next() : 0.10 + 0.15 * uniform.next());
                // Suite 8 Party-Asymmetric: P0 all large, P1 all small
                if (id == 8) weight = (party == 0 ? 3.5 + 1.5 * uniform.next() : 0.20 + 0.35 * uniform.next());
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

        const int id = clampSuiteId(m_suite_id);
        const std::array<size_t, 12> shared_counts = { 2, 2, 1, 2, 1, 2, 2, 2, 20, 6, 8, 16 };
        const size_t desired_shared = shared_counts[static_cast<size_t>(id - 1)];
        auto& uniform = random()->uniform;
        auto unit = [&uniform]() { return uniform.next(); };
        auto tuneForVisibleBasins = [](PartyPeakSpec& peak, int suite_id, size_t dim) {
            // Defaults (used only when a suite doesn't override)
            Real min_radius = 80;
            Real max_condition = 14;
            Real max_rotation = 1.00;
            Real max_linkage = 12;
            Real max_asymmetry = 1.40;
            Real max_deceptive = 0.35;

            // P01 (Suite 1): Smooth balanced
            // Defining feature: MULTIMODAL ENUMERATION in a smooth, near-spherical landscape
            // No rotation/linkage/asymmetry/deceptive pure Gaussian peaks, mild condition
            if (suite_id == 1) {
                min_radius = 80; max_condition = 4.0;
                max_rotation = 0.0; max_linkage = 0.0; max_asymmetry = 0.0; max_deceptive = 0.0;
            }

            // P02 (Suite 2): Unequal k-d tree basin volumes
            // Defining feature: BASIN GEOMETRY â€” some leaves very large, some small
            // Mild per-party condition variation; NO rotation/asymmetry/deceptive
            if (suite_id == 2) {
                min_radius = 55; max_condition = 8.0;
                max_rotation = 0.0; max_linkage = 0.0; max_asymmetry = 0.0; max_deceptive = 0.0;
            }

            // P09 (Suite 3): Neighbor Deception (supplementary) 
            // Narrow shared peaks (h=90) surrounded by broad near-optimal private peaks
            if (suite_id == 3) {
                min_radius = 50; max_condition = 3.0;
                max_rotation = 0.0; max_linkage = 0.0; max_asymmetry = 0.0; max_deceptive = 0.0;
            }

            // P03 (Suite 4): Rugged, asymmetrically transformed subspaces
            // Defining feature: NON-SMOOTH GRADIENT via linkage + asymmetry + deceptive
            // NO rotation (distinguishes P03 from P04/P06)
            if (suite_id == 4) {
                min_radius = 80; max_condition = 3.5;
                max_rotation = 0.0; max_linkage = 12.0; max_asymmetry = 1.40; max_deceptive = 0.32;
            }

            // P04 (Suite 5): Rotated, ill-conditioned basins
            // Defining feature: ORIENTATION + ELONGATION only no ruggedness at all
            if (suite_id == 5) {
                min_radius = 75; max_condition = 16.0; max_rotation = 1.57;
                max_linkage = 0.0; max_asymmetry = 0.0; max_deceptive = 0.0;
            }

            // P10 (Suite 6): Strong Deception (supplementary)
            // Very narrow shared peaks, extremely broad private peaks dominate per-party search
            if (suite_id == 6) {
                min_radius = 40; max_condition = 3.0;
                max_rotation = 0.0; max_linkage = 0.0; max_asymmetry = 0.0; max_deceptive = 0.0;
            }

            // P11 (Suite 7): Multi-scale Hierarchy (supplementary)
            // Extreme leaf-weight contrast; private peaks span two difficulty levels
            if (suite_id == 7) {
                min_radius = 60; max_condition = 22.0; max_rotation = 0.55;
                max_linkage = 0.0; max_asymmetry = 0.0; max_deceptive = 0.10;
            }

            // P12 (Suite 8): Party-Asymmetric (supplementary)
            // P0: broad smooth basins. P1: narrow ill-conditioned rotated basins
            if (suite_id == 8) {
                min_radius = 35; max_condition = 40.0; max_rotation = 0.90;
                max_linkage = 0.0; max_asymmetry = 0.0; max_deceptive = 0.08;
            }

            // P05 (Suite 9): Dense multimodal 20 shared optima
            // Defining feature: ENUMERATION SCALE transforms must be minimal so the
            // challenge is purely locating all 20 optima, not overcoming landscape complexity
            if (suite_id == 9) {
                min_radius = 90; max_condition = 3.5; max_rotation = 0.25;
                max_linkage = 2.0; max_asymmetry = 0.20; max_deceptive = 0.06;
            }

            // P06 (Suite 10): Mixed rotation + ill-cond + rugged + deceptive, 6 opt
            // Strengthened: Îº up to 14, rotation up to 1.2, deceptive up to 0.30
            if (suite_id == 10) {
                min_radius = 70; max_condition = 14.0; max_rotation = 1.20;
                max_linkage = 8.0; max_asymmetry = 0.90; max_deceptive = 0.30;
            }

            // P07 (Suite 11): Hierarchical + disconnected + ill-conditioned, 8 opt
            if (suite_id == 11) {
                min_radius = 65; max_condition = 20.0; max_rotation = 1.10;
                max_linkage = 2.0; max_asymmetry = 0.15; max_deceptive = 0.10;
            }

            // P08 (Suite 12): Dense full-mixed, 16 opt â€” highest overall difficulty
            if (suite_id == 12) {
                min_radius = 60; max_condition = 22.0; max_rotation = 2.00;
                max_linkage = 14.0; max_asymmetry = 1.80; max_deceptive = 0.50;
            }

            // Dimension scaling: Gaussian basin volume (radius/domain_width)^D.
            // Without scaling, a basin that fills 40% of D=2 fills only 0.01% of D=10.
            // Multiplying radius by sqrt(D/2) keeps relative coverage stable across dims.
            const Real dim_scale = std::sqrt(static_cast<Real>(dim) / 2.0);
            peak.radius = std::max(peak.radius * dim_scale, min_radius * dim_scale);
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
                    // P01: smooth near-spherical shared peaks, mild condition only
                    peak.radius = 90 + 40 * unit();
                    peak.condition = 1.0 + 1.5 * unit();  // [1, 2.5]
                }
                if (id == 2) {
                    // P02: party-asymmetric basin widths (P0 wider, P1 narrower, or vice versa)
                    // The leaf WEIGHT already creates volume contrast; radius amplifies the effect
                    if (party == 0) {
                        peak.radius = 120 + 60 * unit();  // broad for party 0
                        peak.condition = 2 + 4 * unit();  // mild condition
                    }
                    else {
                        peak.radius = 55 + 35 * unit();   // narrower for party 1
                        peak.condition = 4 + 6 * unit();  // slightly higher condition
                    }
                }
                if (id == 3) {
                    // Neighbor Deception: shared optima are narrow (small radius) so
                    // the nearby broad private peaks appear more attractive initially
                    peak.radius = 55 + 18 * unit();      // narrow shared basin
                    peak.condition = 1.0 + 0.8 * unit(); // near-spherical â€” clear target
                    // no rotation, no deceptive: challenge is purely spatial
                }
                if (id == 4) {
                    // P03: rugged landscape via linkage + asymmetry + deceptive, NO rotation
                    // Generate within tuneForVisibleBasins caps (max_linkage=12, max_asym=1.4, max_dec=0.32)
                    // so mode diversity survives (previously all values were clipped to identical caps)
                    const size_t mode = (shared_count + 2 * party) % 4;
                    peak.radius = (party == 0 ? 80 + 20 * mode : 90 - 10 * mode) + 20 * unit();
                    peak.linkage = 4.0 + mode * 1.5 + 6.5 * unit();  // [4, 12]
                    peak.asymmetry = 0.50 + mode * 0.15 + 0.75 * unit(); // [0.5, 1.4]
                    peak.deceptive = (mode % 2 == 0 ? 0.08 : 0.16) + 0.16 * unit(); // [0.08, 0.32]
                }
                if (id == 5) {
                    // P04: pure rotation + ill-conditioning NO linkage/asymmetry/deceptive
                    // Parties rotate in opposite directions, creating misaligned elongated ellipsoids
                    peak.radius = 80 + 40 * unit();
                    peak.condition = 6 + 10 * unit();   // [6, 16]
                    peak.rotation = (party == 0 ? 1 : -1) * (0.5 + 1.07 * unit()); // [0.5, 1.57]
                }
                if (id == 6) {
                    // Strong Deception: shared optima have VERY NARROW basins so that
                    // the broad private peaks (h~87) dominate each party's landscape
                    peak.radius = 48 + 12 * unit();       // very narrow hard to locate
                    peak.condition = 1.0 + 0.6 * unit();  // near-spherical, no extra structure
                }
                if (id == 7) {
                    // Multi-scale Hierarchy: shared peaks in SMALL leaves moderate radius,
                    // moderate condition; the challenge is that small leaves are hard to enter
                    peak.radius = 72 + 28 * unit();
                    peak.condition = 3 + 5 * unit();
                }
                if (id == 8) {
                    // Party-Asymmetric: P0 sees a broad, easy basin; P1 sees a narrow, harder one
                    if (party == 0) {
                        peak.radius = 140 + 40 * unit();  // broad, easily approached
                        peak.condition = 2 + 2 * unit();
                    }
                    else {
                        peak.radius = 52 + 18 * unit();   // narrow, harder to locate
                        peak.condition = 10 + 8 * unit(); // moderately ill-conditioned
                    }
                }
                if (id == 9) {
                    // P05: 20 shared optima challenge is ENUMERATION SCALE, not landscape complexity
                    // Generate within caps (max_condition=3.5, max_rotation=0.25, max_linkage=2.0,
                    // max_asymmetry=0.20, max_deceptive=0.06) so mode diversity SURVIVES clipping
                    const size_t mode = shared_count % 5;
                    peak.radius = 85 + 25 * unit();
                    if (mode == 0) {
                        // mild condition + mild opposite-sign rotation
                        peak.condition = 2.0 + 1.5 * unit();  // [2, 3.5]
                        peak.rotation = (party == 0 ? 1 : -1) * (0.08 + 0.17 * unit());  // [0.08, 0.25]
                    }
                    else if (mode == 1) {
                        // mild linkage + asymmetry + deceptive
                        peak.linkage = 0.8 + 1.2 * unit();   // [0.8, 2.0]
                        peak.asymmetry = 0.06 + 0.14 * unit(); // [0.06, 0.20]
                        peak.deceptive = 0.02 + 0.04 * unit(); // [0.02, 0.06]
                    }
                    else if (mode == 2) {
                        // mild condition + mild deceptive gradient
                        peak.condition = 1.8 + 1.7 * unit();   // [1.8, 3.5]
                        peak.deceptive = 0.02 + 0.04 * unit();
                    }
                    else if (mode == 3) {
                        // pure geometry: radius variation, near-spherical
                        peak.radius = 75 + 35 * unit();
                        peak.condition = 1.2 + 1.0 * unit();   // near-spherical
                    }
                    else {
                        // mild rotation + mild deceptive
                        peak.condition = 1.5 + 2.0 * unit();
                        peak.rotation = (party == 0 ? -1 : 1) * (0.08 + 0.17 * unit());
                        peak.deceptive = 0.02 + 0.04 * unit();
                    }
                }
                if (id == 10) {
                    // P06: mixed transforms (max_cond=14, max_rot=1.2, max_link=8,
                    // max_asym=0.90, max_dec=0.30) generate within caps for diversity
                    const size_t mode = shared_count % 6;
                    peak.condition = 5.0 + 9.0 * unit();   // [5, 14]
                    peak.rotation = (party == 0 ? 1 : -1) * (0.4 + 0.8 * unit());  // [0.4, 1.2]
                    peak.radius = (mode % 2 == 0 ? 85 : 60) + 25 * unit();
                    if (mode == 0 || mode == 3) {
                        peak.linkage = 3.0 + 5.0 * unit();   // [3, 8]
                        peak.asymmetry = 0.30 + 0.60 * unit(); // [0.3, 0.9]
                    }
                    if (mode == 1 || mode == 4) {
                        peak.deceptive = 0.10 + 0.20 * unit(); // [0.10, 0.30]
                    }
                    if (mode == 5) {
                        peak.linkage = 2.0 + 5.0 * unit();   // [2, 7]
                        peak.deceptive = 0.08 + 0.14 * unit(); // [0.08, 0.22]
                    }
                }
                if (id == 11) {
                    // P07: ill-conditioned + hierarchical; condition+rotation dominate
                    // (max_cond=20, max_rot=1.10, max_linkage=2, max_asym=0.15, max_dec=0.10)
                    const size_t mode = shared_count % 8;
                    peak.condition = 8.0 + 12.0 * unit();  // [8, 20]
                    peak.rotation = (party == 0 ? -1 : 1) * (0.30 + 0.80 * unit());  // [0.3, 1.1]
                    peak.radius = (mode % 4 == 0 ? 110 : (mode % 4 == 1 ? 55 : 80)) + 20 * unit();
                    // Mild ruggedness for structural variety only
                    if (mode % 4 >= 2) peak.linkage = 0.3 + 1.7 * unit();  // [0.3, 2.0]
                    if (mode % 4 == 3) {
                        peak.asymmetry = 0.04 + 0.11 * unit(); // [0.04, 0.15]
                        peak.deceptive = 0.03 + 0.07 * unit(); // [0.03, 0.10]
                    }
                }
                if (id == 12) {
                    // P08: full-mixed at MAXIMUM levels across all transforms
                    // (max_cond=22, max_rot=2.00, max_linkage=14, max_asym=1.80, max_dec=0.50)
                    const size_t mode = shared_count % 8;
                    peak.condition = 10.0 + 12.0 * unit(); // [10, 22]
                    peak.rotation = (party == 0 ? 1 : -1) * (0.8 + 1.2 * unit());  // [0.8, 2.0]
                    peak.radius = (mode == 0 || mode == 5 ? 90 : (mode == 3 ? 60 : 75)) + 20 * unit();
                    if (mode == 1 || mode == 3 || mode == 7) {
                        peak.deceptive = 0.18 + 0.32 * unit();  // [0.18, 0.50]
                    }
                    if (mode == 2 || mode == 4 || mode == 7) {
                        peak.linkage = 6.0 + 8.0 * unit();    // [6, 14]
                        peak.asymmetry = 0.70 + 1.10 * unit();  // [0.70, 1.80]
                    }
                }
                tuneForVisibleBasins(peak, id, m_problem_dimension);
                peak.height = Real(90);  // shared optima always at exactly 90
                spec.parties[party].peaks.push_back(peak);
            }
            ++shared_count;
        }

        // (p6_neighbor_deceptive_done removed suite 6 redesigned as Strong Deception)
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
                    // P01 private peaks: smooth near-spherical, mild condition only
                    // Private peaks must be attractive per-party (h 70-85) but Q=min_p F_p
                    // is low at these positions mediating population gravitates to shared optima
                    peak.height = 65 + 22 * unit();
                    peak.radius = 90 + 50 * unit();    // moderate radius, always min_radius=80
                    peak.condition = 1.0 + 1.8 * unit(); // [1, 2.8] very mild
                    // rotation=0, linkage=0, asymmetry=0, deceptive=0 (smooth landscape)
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
                    // Neighbor Deception: private peaks are TALL and BROAD, centered near
                    // the shared optimum â€” creating a "shadow attractor"
                    peak.height = 82 + 6 * unit();     // near-optimal h [82, 88)
                    peak.radius = 220 + 60 * unit();   // very broad basin, overlaps shared region
                    peak.condition = 1.0 + 1.0 * unit();
                    // no rotation, no deceptive: spatial proximity IS the deception
                }
                if (id == 5) {
                    // P04 private peaks: rotation + condition only matches strengthened caps
                    peak.radius = 80 + 50 * unit();
                    peak.condition = 5 + 11 * unit();              // [5, 16]
                    peak.rotation = (party == 0 ? 1 : -1) * (0.4 + 1.17 * unit()); // [0.4, 1.57]
                    // linkage/asymmetry/deceptive intentionally left at 0
                }
                if (id == 4) {
                    // P03 private peaks: same rugged profile as shared generate within caps
                    const size_t mode = (i + 3 * party) % 5;
                    peak.radius = 80 + 40 * unit();
                    peak.linkage = 4.0 + 8.0 * unit();              // [4, 12]
                    peak.asymmetry = 0.50 + 0.90 * unit();            // [0.5, 1.4]
                    peak.deceptive = (mode % 2 == 0 ? 0.08 : 0.16) + 0.16 * unit(); // [0.08, 0.32]
                }
                if (id == 6) {
                    // Strong Deception: private peaks are very broad and near-optimal (h~87)
                    // so that single-party search finds THESE, not the narrow h=90 shared peaks
                    peak.height = 84 + 4 * unit();      // h [84, 88) near-optimal
                    peak.radius = 260 + 80 * unit();    // extremely broad global attractor
                    peak.condition = 1.0 + 0.6 * unit(); // spherical easy to climb
                }
                if (id == 7) {
                    // Multi-scale Hierarchy: alternate between very broad (Level-1) peaks
                    // that guide search into a region, and narrower (Level-2) peaks inside
                    const size_t level = i % 2;
                    if (level == 0) {
                        // Level-1: very broad outer basin draws search to the region
                        peak.height = 62 + 16 * unit();
                        peak.radius = 250 + 60 * unit();
                        peak.condition = 2 + 2 * unit();
                    }
                    else {
                        // Level-2: narrower inner basin funnel toward shared optimum vicinity
                        peak.height = 52 + 24 * unit();
                        peak.radius = 80 + 40 * unit();
                        peak.condition = 8 + 12 * unit();
                    }
                }
                if (id == 8) {
                    // Party-Asymmetric: party-0 landscape is broad and easy;
                    // party-1 landscape is narrow, ill-conditioned and rotated
                    if (party == 0) {
                        peak.height = 60 + 20 * unit();
                        peak.radius = 190 + 50 * unit();  // broad, easy
                        peak.condition = 2 + 3 * unit();
                    }
                    else {
                        peak.height = 54 + 24 * unit();
                        peak.radius = 40 + 28 * unit();   // narrow, hard to find
                        peak.condition = 18 + 20 * unit(); // strongly ill-conditioned
                        peak.rotation = (i % 2 == 0 ? 1 : -1) * (0.7 + 0.8 * unit());
                    }
                }
                if (id == 10) {
                    // P06 private peaks: mixed transforms within caps (same as shared)
                    const size_t mode = (i * 3 + party) % 7;
                    peak.height = party == 0 ? 60 + 24 * unit() : 54 + 32 * unit();
                    peak.condition = 5.0 + 9.0 * unit();   // [5, 14]
                    peak.rotation = (party == 0 ? 1 : -1) * (0.4 + 0.8 * unit());  // [0.4, 1.2]
                    peak.radius = (mode == 0 || mode == 4 ? 95 : (mode == 2 ? 60 : 75)) + 20 * unit();
                    if (mode == 2 || mode == 6) {
                        peak.deceptive = 0.10 + 0.20 * unit(); // [0.10, 0.30]
                    }
                    if (mode == 0 || mode == 3 || mode == 5) {
                        peak.linkage = 3.0 + 5.0 * unit();   // [3, 8]
                        peak.asymmetry = 0.30 + 0.60 * unit(); // [0.3, 0.9]
                    }
                }
                if (id == 11) {
                    // P07 private peaks: ill-conditioned + hierarchical within caps
                    // (max_cond=20, max_rot=1.10, max_linkage=2.0, max_asym=0.15, max_dec=0.10)
                    const size_t mode = (i + 2 * party) % 8;
                    peak.height = party == 0 ? 57 + 27 * unit() : 62 + 22 * unit();
                    peak.condition = 9.0 + 11.0 * unit();  // [9, 20]
                    peak.rotation = (party == 0 ? -1 : 1) * (0.30 + 0.80 * unit());  // [0.3, 1.1]
                    peak.radius = (mode < 2 ? 110 : (mode == 3 || mode == 6 ? 60 : 80)) + 20 * unit();
                    if (mode % 4 >= 2) peak.linkage = 0.3 + 1.7 * unit();  // [0.3, 2.0]
                    if (mode % 4 == 3) {
                        peak.asymmetry = 0.04 + 0.11 * unit(); // [0.04, 0.15]
                        peak.deceptive = 0.03 + 0.07 * unit(); // [0.03, 0.10]
                    }
                }
                if (id == 12) {
                    // P08 private peaks: full-mixed maximum transforms within caps
                    // (max_cond=22, max_rot=2.00, max_linkage=14, max_asym=1.80, max_dec=0.50)
                    const size_t mode = (i * 5 + party * 3) % 10;
                    peak.height = party == 0 ? 54 + 34 * unit() : 50 + 38 * unit();
                    peak.condition = 10.0 + 12.0 * unit(); // [10, 22]
                    peak.rotation = (party == 0 ? 1 : -1) * (0.7 + 1.3 * unit());  // [0.7, 2.0]
                    peak.radius = (mode == 0 || mode == 7 ? 100 : (mode == 2 || mode == 5 ? 60 : 75)) + 20 * unit();
                    if (mode == 3 || mode == 5 || mode == 9) {
                        peak.deceptive = 0.18 + 0.32 * unit();  // [0.18, 0.50]
                    }
                    if (mode == 0 || mode == 4 || mode == 6 || mode == 8) {
                        peak.linkage = 6.0 + 8.0 * unit();    // [6, 14]
                        peak.asymmetry = 0.70 + 1.10 * unit();  // [0.70, 1.80]
                    }
                }
                tuneForVisibleBasins(peak, id, m_problem_dimension);
                if (peak.height >= Real(90)) peak.height = Real(89);  // private peaks must not reach 90
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