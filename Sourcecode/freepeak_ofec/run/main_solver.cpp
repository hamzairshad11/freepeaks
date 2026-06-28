#include "mpmcoea_solver.hpp"
#include "custom_method.hpp"
#include "instance/problem/continuous/free_peaks/factory.h"

int main() {
    try {
        // Step 1: Run the algorithm experiment (all 12 problems, 4 dims)
        std::cout << ">>> [MPM  M-CoEA] Starting Experiment with Visualization\n" << std::endl;
        ofec::runMPMMCoEAExperiment();

        // Step 2: Generate high-quality benchmark landscape files (D=2 only)
        // Produces grid_2d.tsv (400x400) + optima_2d.tsv for each of the 12 suites.
        std::cout << "\n>>> [Landscape] Generating 400x400 benchmark landscapes for P1-P12...\n";
        const std::string thesis_root =
            "E:\\HITSZ\\Research\\Multimodal_Multiparty_Optimization\\ThesisProject\\";
        outputFreePeaksMultipartyLandscape12(thesis_root);
        std::cout << ">>> [Landscape] Done.\n";
    }
    catch (const std::exception& e) {
        std::cerr << "[ERROR] " << e.what() << std::endl;
        return 1;
    }
    return 0;
}