#ifndef OFEC_FREE_PEAKS_CONSTRAINT_BASE_H
#define OFEC_FREE_PEAKS_CONSTRAINT_BASE_H

#include "../../../../../../core/parameter/parameter_map.h"
#include "../distance/distance_base.h"
#include "../../type.h"


namespace ofec {
	namespace free_peaks {
		class ConstraintBase {
		protected:
			Problem* m_pro = nullptr;
			const std::string m_subspace_name;
			const ParameterMap m_param;
			const DistanceBase* m_dis;

		public:
			ConstraintBase(Problem *pro, const std::string& subspace_name, const ParameterMap& param);
			virtual ~ConstraintBase() = default;
			virtual double evaluate(const std::vector<Real>& var) const = 0;
			virtual void bindData();
			const ParameterMap& parameter() const { return m_param; }
		};
	}
}

#endif // !OFEC_FREE_PEAKS_CONSTRAINT_BASE_H
