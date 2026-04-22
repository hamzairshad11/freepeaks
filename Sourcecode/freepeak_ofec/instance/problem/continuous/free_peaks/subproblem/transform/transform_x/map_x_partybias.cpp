#include "map_x_partybias.h"
#include "../../../free_peaks.h" 

#include <cmath>
#include <vector>

namespace ofec {
    namespace free_peaks {

        //  addInputParameters
        //  Called once by the macro-generated constructor via:
        //    OFEC_ABSTRACT_INSTANCE  ->  constructor  ->  addInputParameters().
        //  X_TransformBase's own macro-generated constructor already
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

        //  initialize
        //  Framework pattern followed by every Transform subclass:
        //    1. Call TransformBase::initialize()
        //       -> sets m_pro, m_subspace_name
        //       -> calls m_input_parameters.input(param), which copies every
        //          value from the ParameterMap into the registered member refs
        //          (m_party_id, m_magnitude, m_rotation_angle, m_condition_number).
        //    Use those now-populated members to build the derived matrices.
        void MapXPartyBias::initialize(Problem* pro,
            const std::string& subspace_name,
            const ParameterMap& param)
        {
            // Base initialisation (also populates member vars)
            TransformBase::initialize(pro, subspace_name, param);

            const int dim = static_cast<int>(pro->numberVariables());

            // Party-specific bias (shift) vector
            // Deterministic and unique per (party_id, dimension) pair.
            // A sin phase difference avoids axis-aligned symmetry between parties.
            m_shift.assign(static_cast<size_t>(dim), Real(0));
            for (int i = 0; i < dim; ++i) {
                m_shift[static_cast<size_t>(i)] =
                    m_magnitude * std::sin(static_cast<Real>(m_party_id) * Real(2.1)
                        + static_cast<Real>(i) * Real(1.3));
            }

            // Rotation matrix (Givens rotation in the 0-1 plane) --
            // For dim > 2 the remaining axes are kept at identity.
            // Full random orthogonal matrix
            m_rotation_matrix.randomize(&pro->random()->uniform);
            m_rotation_matrix.generateRotationClassical(&pro->random()->normal, 1.0);
            m_rotation_matrix.resize(static_cast<size_t>(dim),
                static_cast<size_t>(dim));
            m_rotation_matrix.identify();

            if (dim >= 2 && m_rotation_angle != Real(0)) {
                const Real c = std::cos(m_rotation_angle);
                const Real s = std::sin(m_rotation_angle);
                m_rotation_matrix[0][0] = c;
                m_rotation_matrix[0][1] = -s;
                m_rotation_matrix[1][0] = s;
                m_rotation_matrix[1][1] = c;
            }

            // Diagonal ill-conditioning matrix
            // Scales axis i by kappa^(i/(d-1)), giving an effective condition
            // number of kappa == m_condition_number.
            // BBOB-style ill-conditioning convention
            m_scaling_matrix.resize(static_cast<size_t>(dim),
                static_cast<size_t>(dim));
            m_scaling_matrix.identify();

            if (m_condition_number > Real(1) && dim >= 2) {
                const Real exponent = Real(1) / static_cast<Real>(dim - 1);
                const Real base = std::pow(m_condition_number, exponent);
                for (int i = 0; i < dim; ++i) {
                    m_scaling_matrix[static_cast<size_t>(i)][static_cast<size_t>(i)] =
                        std::pow(base, static_cast<Real>(i));
                }
            }

            m_initialized = true;
        }

        //  transfer  (core interface)
        //  Order: bias -> rotation -> ill-conditioning.
        //  Each step pivots around the party's own shift point so the operations
        //  compose geometrically.
        void MapXPartyBias::transfer(std::vector<Real>& x,
            const std::vector<Real>& /*var*/)
        {
            if (!m_initialized || x.empty()) return;

            applyBias(x);
            applyRotation(x);
            applyIllConditioning(x);
        }

        //  applyBias
        void MapXPartyBias::applyBias(std::vector<Real>& x) const {
            for (size_t i = 0; i < x.size(); ++i) {
                x[i] += m_shift[i];
                clampToUnitInterval(x[i]);
            }
        }

        //  applyRotation
        //  y = R * (x - shift) + shift
        void MapXPartyBias::applyRotation(std::vector<Real>& x) const {
            const int n = static_cast<int>(x.size());

            std::vector<Real> centred(static_cast<size_t>(n));
            for (int i = 0; i < n; ++i)
                centred[static_cast<size_t>(i)] = x[static_cast<size_t>(i)]
                - m_shift[static_cast<size_t>(i)];

            std::vector<Real> rotated(static_cast<size_t>(n), Real(0));
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    rotated[static_cast<size_t>(i)] +=
                    m_rotation_matrix[static_cast<size_t>(i)][static_cast<size_t>(j)]
                    * centred[static_cast<size_t>(j)];

            for (int i = 0; i < n; ++i) {
                rotated[static_cast<size_t>(i)] += m_shift[static_cast<size_t>(i)];
                clampToUnitInterval(rotated[static_cast<size_t>(i)]);
            }

            x = std::move(rotated);
        }

        //  applyIllConditioning
        //  y_i = S_ii * (x_i - shift_i) + shift_i   (O(n) because S is diagonal)
        void MapXPartyBias::applyIllConditioning(std::vector<Real>& x) const {
            const int n = static_cast<int>(x.size());

            std::vector<Real> scaled(static_cast<size_t>(n));
            for (int i = 0; i < n; ++i) {
                const size_t ui = static_cast<size_t>(i);
                scaled[ui] = m_scaling_matrix[ui][ui] * (x[ui] - m_shift[ui])
                    + m_shift[ui];
                clampToUnitInterval(scaled[ui]);
            }

            x = std::move(scaled);
        }

        //  clampToUnitInterval  (static helper)
        //  Clamp rather than wrap: wrapping creates C0 discontinuities that
        //  distort landscape topology and mislead heuristic search.
        /*static*/ void MapXPartyBias::clampToUnitInterval(Real& val) {
            if (val < Real(0)) val = Real(0);
            else if (val > Real(1)) val = Real(1);
        }

    } // namespace free_peaks
} // namespace ofec