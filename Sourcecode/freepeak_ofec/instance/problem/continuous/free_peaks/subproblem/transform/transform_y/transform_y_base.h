#ifndef OFEC_FREE_PEAKS_TRANSFORM_BASE_Y_H
#define OFEC_FREE_PEAKS_TRANSFORM_BASE_Y_H

#include "../transform_base.h"

namespace ofec {
	class Problem;
	namespace free_peaks {
		class Subproblem;
		class  Y_TransformBase : public TransformBase {
			OFEC_ABSTRACT_INSTANCE(Y_TransformBase)
		private:
			std::string m_register_name = "Y_TransformBase";
		public:

			void setRegisterName(const std::string& name) {
				m_register_name = name;
			}

			void addInputParameters() {
				m_input_parameters.add("transformY", new InputParameterValueTypeString(m_register_name));
			}
		};
	}
}

#endif // !OFEC_FREE_PEAKS_TRANSFORM_BASE_H