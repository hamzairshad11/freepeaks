#include "map_x_partybias.h"
#include "instance/problem/continuous/free_peaks/free_peaks.h"

namespace ofec {
    namespace free_peaks {

        MapXPartyBias::MapXPartyBias() = default;
        MapXPartyBias::~MapXPartyBias() = default;

        void MapXPartyBias::initialize(Problem* prob, const std::string& subspace_name, const ParameterMap& param) {
            // Call base initialize
            X_TransformBase::initialize(prob, subspace_name, param);

            // Extract party-specific parameters using has() and get()
            if (param.has("party_id")) {
                m_party_id = param.get<int>("party_id");
            }

            if (param.has("magnitude")) {
                m_magnitude = param.get<Real>("magnitude");
            }
        }

        void MapXPartyBias::transform(VariableVector<Real>& x, const std::string& subspace_name) {
            // Apply party-specific bias with wrapping to [0,1]
            auto& vec = x.vector();
            // Create party-specific offset (different for each party)
            Real party_offset = static_cast<Real>(m_party_id) * m_magnitude;

            // Add index variable i
            for (size_t i = 0; i < vec.size(); ++i) {
                auto& val = vec[i];
                Real dim_bias = party_offset * (i + 1);  // Different per dimension
                val += dim_bias;  // Fixed: was adding m_magnitude, should add dim_bias
                // Wrap to [0, 1] range
                while (val < static_cast<Real>(0.0)) val += static_cast<Real>(1.0);
                while (val > static_cast<Real>(1.0)) val -= static_cast<Real>(1.0);
            }
        }

        std::string MapXPartyBias::getName() const {
            return "MapXPartyBias";
        }

    } // namespace free_peaks
} // namespace ofec