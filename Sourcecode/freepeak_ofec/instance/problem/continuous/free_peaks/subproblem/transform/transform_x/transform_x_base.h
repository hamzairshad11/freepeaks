#ifndef OFEC_FREE_PEAKS_TRANSFORM_BASE_X_H
#define OFEC_FREE_PEAKS_TRANSFORM_BASE_X_H

#include "../transform_base.h"

namespace ofec {
	class Problem;
	namespace free_peaks {
		class Subproblem;
		class  X_TransformBase : public TransformBase {
			OFEC_ABSTRACT_INSTANCE(X_TransformBase)
		private:
			std::string m_register_name = "X_TransformBase";
		public:

			void setRegisterName(const std::string& name) {
				m_register_name = name;
			}

			void addInputParameters() {
				m_input_parameters.add("transformX", new InputParameterValueTypeString(m_register_name));
			}
		};
	}
}

#endif // !OFEC_FREE_PEAKS_TRANSFORM_BASE_H