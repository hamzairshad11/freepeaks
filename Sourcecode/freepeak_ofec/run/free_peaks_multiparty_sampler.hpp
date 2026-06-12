#ifndef OFEC_FREE_PEAKS_MULTIPARTY_SAMPLER_HPP
#define OFEC_FREE_PEAKS_MULTIPARTY_SAMPLER_HPP

#include "interface.h"
#include "../core/global.h"
#include "../core/problem/solution.h"
#include "../instance/problem/continuous/free_peaks/free_peaks_multiparty.h"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

namespace ofec {

    inline void sampleFreePeaksMultiPartySuites(
        const std::string& output_dir = "Visualization/free_peaks_multiparty_samples",
        int resolution = 241)
    {
        registerInstance();

        if (resolution < 2) resolution = 2;
        std::filesystem::create_directories(output_dir);

        for (int suite_id = 1; suite_id <= 12; ++suite_id) {
            std::shared_ptr<Environment> env(generateEnvironmentByFactory("free_peaks_multiparty"));
            env->setProblem(generateProblemByFactory("free_peaks_multiparty"));

            ParameterMap pm;
            pm["suite_id"] = suite_id;
            env->problem()->inputParameters().input(pm);
            env->initializeProblem(0.5);

            auto* problem = dynamic_cast<FreePeaksMultiParty*>(env->problem());
            if (!problem) throw std::runtime_error("free_peaks_multiparty initialization failed");

            const auto& spec = problem->currentSpec();
            const std::filesystem::path suite_dir = std::filesystem::path(output_dir) /
                ("suite_" + std::to_string(suite_id) + "_" + spec.name);
            std::filesystem::create_directories(suite_dir);

            std::ofstream grid_out(suite_dir / "grid.txt");
            grid_out << std::fixed << std::setprecision(10);
            grid_out << "# suite_id " << suite_id << "\n";
            grid_out << "# name " << spec.name << "\n";
            grid_out << "# feature " << spec.feature << "\n";
            grid_out << "# columns x0_norm x1_norm x0_plot x1_plot f_party0 f_party1 consensus_min\n";

            auto sol = std::unique_ptr<SolutionBase>(env->problem()->createSolution());
            auto& x = dynamic_cast<Solution<VariableVector<Real>>&>(*sol).variable().vector();
            for (int iy = 0; iy < resolution; ++iy) {
                const Real x1 = static_cast<Real>(iy) / static_cast<Real>(resolution - 1);
                for (int ix = 0; ix < resolution; ++ix) {
                    const Real x0 = static_cast<Real>(ix) / static_cast<Real>(resolution - 1);
                    x[0] = x0;
                    x[1] = x1;
                    sol->evaluate(env.get(), false);
                    const Real f0 = sol->objective(0);
                    const Real f1 = sol->objective(1);
                    grid_out << x0 << ' ' << x1 << ' '
                             << (x0 * 7.0 - 3.5) << ' ' << (x1 * 7.0 - 3.5) << ' '
                             << f0 << ' ' << f1 << ' ' << std::min(f0, f1) << "\n";
                }
                grid_out << "\n";
            }

            std::ofstream opt_out(suite_dir / "shared_optima.txt");
            opt_out << std::fixed << std::setprecision(10);
            opt_out << "# columns x0_norm x1_norm x0_plot x1_plot f_party0 f_party1 consensus_min\n";
            std::ofstream region_out(suite_dir / "shared_regions.txt");
            region_out << "# columns x0_plot x1_plot party0_region party1_region\n";
            for (const auto& pt : spec.shared_optima) {
                x[0] = pt[0];
                x[1] = pt[1];
                sol->evaluate(env.get(), false);
                const Real f0 = sol->objective(0);
                const Real f1 = sol->objective(1);
                opt_out << pt[0] << ' ' << pt[1] << ' '
                        << (pt[0] * 7.0 - 3.5) << ' ' << (pt[1] * 7.0 - 3.5) << ' '
                        << f0 << ' ' << f1 << ' ' << std::min(f0, f1) << "\n";
                std::vector<Real> pt_vec = { pt[0], pt[1] };
                region_out << (pt[0] * 7.0 - 3.5) << ' ' << (pt[1] * 7.0 - 3.5) << ' '
                           << problem->partyRegionName(0, pt_vec) << ' '
                           << problem->partyRegionName(1, pt_vec) << "\n";
            }

            std::cout << "[sample] suite " << suite_id << " -> "
                      << suite_dir.string() << "\n";
        }
    }
}

#endif
