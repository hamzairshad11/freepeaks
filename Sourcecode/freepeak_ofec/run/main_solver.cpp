#include "mpmcoea_solver.hpp"
#include "instance/problem/continuous/free_peaks/factory.h"
#include <iostream>

int main() {
    try {
        ofec::runMPMCoEAExperiment();
    }
    catch (const std::exception& e) {
        std::cerr << "[ERROR] " << e.what() << std::endl;
        return 1;
    }
    return 0;
}