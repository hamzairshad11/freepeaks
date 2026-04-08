#include "mpmcoea_solver.hpp"
#include "../instance/problem/continuous/free_peaks/subproblem/transform/transform_x/map_x_partybias.h"
#include "instance/problem/continuous/free_peaks/factory.h"
#include <iostream>

namespace ofec {
    namespace free_peaks {

        void registerCustomTransforms() {
            // Register the MapXPartyBias transform with the factory
            FactoryFP<X_TransformBase>::REGISTER(
                "MapXPartyBias",
                []() { return std::make_unique<MapXPartyBias>(); }
            );

            std::cout << "[INFO] Custom transforms registered successfully" << std::endl;
        }

    } // namespace free_peaks
} // namespace ofec

int main() {
    // Register custom transforms before running experiment
    ofec::free_peaks::registerCustomTransforms();

    try {
        ofec::runMPMCoEAExperiment();
    }
    catch (const std::exception& e) {
        std::cerr << "[ERROR] " << e.what() << std::endl;
        return 1;
    }
    return 0;
}