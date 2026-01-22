#ifndef OFEC_FREE_PEAKS_TRANSFORM_BASE_H
#define OFEC_FREE_PEAKS_TRANSFORM_BASE_H

#include "../../../../../../core/parameter/parameter_map.h"
#include "../../type.h"

namespace ofec {
	class Problem;
	namespace free_peaks {
		class Subproblem;
		class  TransformBase {
		protected:
			Problem* m_pro = nullptr;
			const std::string m_subspace_name;
			const ParameterMap m_param;
			Subproblem* m_sub_pro = nullptr;
		public:
			TransformBase(Problem *pro, const std::string& subspace_name, const ParameterMap& param);
			virtual ~TransformBase() = default;
			virtual void transfer(std::vector<Real>& vec, const std::vector<Real>& var) = 0;
			virtual void bindData() = 0;
			virtual void bindSubPro(Subproblem* subpro) {
				m_sub_pro = subpro;
			}
			const ParameterMap& parameter() const { return m_param; }
		};
	}
}

#endif // !OFEC_FREE_PEAKS_TRANSFORM_BASE_H