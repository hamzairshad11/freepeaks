#include "map_x_partybias.h"
#include "../../../free_peaks.h" 

#include <cmath>
#include <vector>

namespace ofec {
    namespace free_peaks {

        
        //OFEC_REGISTER_INSTANCE(MapXPartyBias, ofec::free_peaks::MapXPartyBias)

        void MapXPartyBias::addInputParameters() {
            m_input_parameters.add(
                "party_id",
                new InputParameterValueTypeInt(m_party_id));

            m_input_parameters.add(
                "magnitude",
                new InputParameterValueTypeReal(m_magnitude));

            m_input_parameters.add(
                "rotation_angle",
                new InputParameterValueTypeReal(m_rotation_angle));

            m_input_parameters.add(
                "condition_number",
                new InputParameterValueTypeReal(m_condition_number));
        }

        void MapXPartyBias::initialize(Problem* pro,
            const std::string& subspace_name,
            const ParameterMap& param)
        {
            // Base initialisation (populates m_party_id, m_magnitude, etc.)
            TransformBase::initialize(pro, subspace_name, param);

            const int dim = static_cast<int>(pro->numberVariables());

            // Bias vector: deterministic per (party_id, dimension)
            m_shift.assign(dim, Real(0));
            for (int i = 0; i < dim; ++i) {
                double phase = static_cast<double>(m_party_id) * 2.1
                    + static_cast<double>(i) * 1.3;
                m_shift[i] = m_magnitude * static_cast<Real>(std::sin(phase));
            }

            // Rotation matrix
            // Start with identity. If rotation_angle != 0, apply a Givens
            // rotation in the plane (party_id % dim, (party_id+1) % dim).
            // This ensures Party 0 and Party 1 couple *different* axes when dim>=3.
            m_rotation_matrix.resize(dim, dim);
            m_rotation_matrix.identify();

            if (dim >= 2 && std::abs(m_rotation_angle) > Real(1e-12)) {
                int d1 = m_party_id % dim;
                int d2 = (m_party_id + 1) % dim;
                Real c = std::cos(m_rotation_angle);
                Real s = std::sin(m_rotation_angle);

                m_rotation_matrix[d1][d1] = c;
                m_rotation_matrix[d1][d2] = -s;
                m_rotation_matrix[d2][d1] = s;
                m_rotation_matrix[d2][d2] = c;
            }

            // Ill-conditioning (diagonal scaling)
            // BBOB convention: axis i scaled by kappa^(i/(d-1)).
            // If condition_number <= 1, we leave it as identity (no stretch).
            m_scaling_matrix.resize(dim, dim);
            m_scaling_matrix.identify();

            if (m_condition_number > Real(1) && dim >= 2) {
                Real exponent = Real(1) / static_cast<Real>(dim - 1);
                Real base = std::pow(m_condition_number, exponent);
                for (int i = 0; i < dim; ++i) {
                    m_scaling_matrix[i][i] = std::pow(base, static_cast<Real>(i));
                }
            }

            m_initialized = true;
        }

        // Core interface: maps x in-place.
        // Net mathematical effect: x_out = S * R * x_in + shift
        void MapXPartyBias::transfer(std::vector<Real>& x,
            const std::vector<Real>& /*var*/)
        {
            if (!m_initialized || x.empty()) return;
            applyBias(x);            // x = x_in + shift
            applyRotation(x);        // x = R*(x_in) + shift
            applyIllConditioning(x); // x = S*R*x_in + shift
        }

        void MapXPartyBias::applyBias(std::vector<Real>& x) const {
            for (size_t i = 0; i < x.size(); ++i) x[i] += m_shift[i];
        }

        void MapXPartyBias::applyRotation(std::vector<Real>& x) const {
            const int n = static_cast<int>(x.size());
            // Because applyBias was already called, x currently equals x_original + shift.
            // We pivot around the shift point: y = R*(x - shift) + shift.
            std::vector<Real> centred(n);
            for (int i = 0; i < n; ++i) centred[i] = x[i] - m_shift[i];

            std::vector<Real> rotated(n, Real(0));
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    rotated[i] += m_rotation_matrix[i][j] * centred[j];

            for (int i = 0; i < n; ++i) x[i] = rotated[i] + m_shift[i];
        }

        void MapXPartyBias::applyIllConditioning(std::vector<Real>& x) const {
            const int n = static_cast<int>(x.size());
            // Pivot around shift: y = S*(x - shift) + shift
            for (int i = 0; i < n; ++i) {
                x[i] = m_scaling_matrix[i][i] * (x[i] - m_shift[i]) + m_shift[i];
            }
        }

    } // namespace free_peaks
} // namespace ofec