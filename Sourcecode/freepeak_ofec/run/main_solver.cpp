#include "mpmcoea_solver.hpp"
#include "free_peaks_multiparty_sampler.hpp"
#include "../instance/problem/continuous/free_peaks/factory.h"

int main(int argc, char* argv[]) {
    try {
        if (argc > 1 && std::string(argv[1]) == "sample-multiparty") {
            ofec::sampleFreePeaksMultiPartySuites();
            return 0;
        }
        std::cout << ">>> [MPM-CoEA] Starting Experiment\n";
        ofec::runMPMCoEAExperiment();
    }
    catch (const std::exception& e) {
        std::cerr << "[ERROR] " << e.what() << std::endl;
        return 1;
    }
    return 0;
}