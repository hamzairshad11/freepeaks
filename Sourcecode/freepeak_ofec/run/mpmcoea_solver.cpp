#include "mpmcoea_solver.hpp"

namespace ofec {
    namespace free_peaks {

        void registerCustomTransforms() {
            // MapXPartyBias is already registered via REGISTER_FP macro
            // in register_transform_x.cpp, so no manual registration needed here.
            // This function is kept for compatibility but does nothing.
        }

    } // namespace free_peaks
} // namespace ofec