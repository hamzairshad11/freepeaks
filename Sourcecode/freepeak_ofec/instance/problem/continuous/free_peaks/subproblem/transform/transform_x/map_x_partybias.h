//#ifndef OFEC_FREE_PEAKS_MAP_X_PARTY_BIAS_H
//#define OFEC_FREE_PEAKS_MAP_X_PARTY_BIAS_H

#include "transform_x_base.h"
#include "../transform_base.h" 
#include "utility/linear_algebra/matrix.h"
#define _USE_MATH_DEFINES
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
namespace ofec {
    namespace free_peaks {

        class MapXPartyBias : public X_TransformBase {
            OFEC_CONCRETE_INSTANCE(MapXPartyBias)

        protected:
            int m_party_id;
            double m_magnitude;
            double m_rotation_angle;
            double m_condition_number;

            std::vector<Real> m_shift;
            Matrix m_rotation_matrix;
            Matrix m_scaling_matrix;
            bool m_initialized;

        public:
            MapXPartyBias() : m_party_id(0), m_magnitude(0.0), m_rotation_angle(0.0),
                m_condition_number(1.0), m_initialized(false) {
            }
            
            void addInputParameters() {
                X_TransformBase::addInputParameters();
                m_input_parameters.add("party_id", new InputParameterValueTypeInt(m_party_id));
                m_input_parameters.add("magnitude", new InputParameterValueTypeReal(m_magnitude));
                m_input_parameters.add("rotation_angle", new InputParameterValueTypeReal(m_rotation_angle));
                m_input_parameters.add("condition_number", new InputParameterValueTypeReal(m_condition_number));
            }


            void initialize(Problem* prob, const std::string& subspace_name, const ParameterMap& param) override;
            void transfer(std::vector<Real>& x, const std::vector<Real>& var);


        private:
            void applyBias(std::vector<Real>& x) const;
            void applyRotation(std::vector<Real>& x) const;
            void applyIllConditioning(std::vector<Real>& x) const;
            void wrapToUnitInterval(Real& val) const;
        };

    } // namespace free_peaks
} // namespace ofec

//#endif // OFEC_FREE_PEAKS_MAP_X_PARTY_BIAS_H
