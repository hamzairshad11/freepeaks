#include "mpmcoea_solver.hpp"
#include "instance/problem/continuous/free_peaks/factory.h" 

int main() {
    try {
        std::cout << ">>> [MPM-CoEA] Starting Experiment with Visualization\n" << std::endl;

        // Run the experiment
        ofec::runMPMCoEAExperiment();
    }
    catch (const std::exception& e) {
        std::cerr << "[ERROR] " << e.what() << std::endl;
        return 1;
    }
    return 0;
}