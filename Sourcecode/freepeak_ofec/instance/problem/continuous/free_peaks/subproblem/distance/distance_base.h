#ifndef OFEC_FREE_PEAKS_DISTANCE_BASE_H
#define OFEC_FREE_PEAKS_DISTANCE_BASE_H

#include "../../../../../../core/parameter/parameter_map.h"
#include <vector>
#include "../../../../../../core/environment/environment.h"

namespace ofec {
	class Problem;
	namespace free_peaks {

		class  DistanceBase {
		protected:
			Problem* m_pro = nullptr;
			const std::string m_subspace_name;
			const ParameterMap m_param;
			static const Real ms_infinity;

		public:
			DistanceBase(Problem *pro, const std::string& subspace_name, const ParameterMap& param);
			virtual ~DistanceBase() = default;
			virtual Real operator()(const std::vector<Real>& a, const std::vector<Real>& b) const = 0;
			virtual void bindData() {}
			const ParameterMap& parameter() const { return m_param; }
			Real infiniy() const { return ms_infinity; }
		};
	}
}

#endif // !OFEC_FREE_PEAKS_DISTANCE_BASE_H