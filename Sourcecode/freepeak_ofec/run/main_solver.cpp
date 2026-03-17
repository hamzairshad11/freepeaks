#include "mpmcoea_solver.hpp"

int main() {
    try {
        ofec::runMPMCoEAExperiment();  // Runs solver on PRE-GENERATED benchmark
    }
    catch (const std::exception& e) {
        std::cerr << "[ERROR] " << e.what() << std::endl;
        return 1;
    }
    return 0;
}