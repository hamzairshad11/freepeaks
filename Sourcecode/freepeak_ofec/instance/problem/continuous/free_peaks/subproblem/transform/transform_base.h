#ifndef OFEC_FREE_PEAKS_TRANSFORM_BASE_H
#define OFEC_FREE_PEAKS_TRANSFORM_BASE_H

#include "../../../../../../core/parameter/parameter_map.h"
#include "../../../../../../core/instance.h"
#include "../../type.h"

namespace ofec {
	class Problem;
	namespace free_peaks {
		class Subproblem;
		class  TransformBase : virtual public Instance {
			OFEC_ABSTRACT_INSTANCE(TransformBase)

		protected:


			Problem* m_pro = nullptr;
			std::string m_subspace_name;
			Subproblem* m_sub_pro = nullptr;
		public:

			virtual void bindData();

			void addInputParameters() {
				//m_input_parameters.add("transform", new InputParameterValueTypeString(m_register_name));
			}	
			virtual void initialize(Problem *pro, const std::string& subspace_name, const ParameterMap& param);
			virtual void transfer(std::vector<Real>& vec, const std::vector<Real>& var) = 0;
			//const ParameterMap& parameter() const { return m_param; }
		};
	}
}

#endif // !OFEC_FREE_PEAKS_TRANSFORM_BASE_H