#ifndef OFEC_FREE_PEAKS_LINEAR_MAP_H
#define OFEC_FREE_PEAKS_LINEAR_MAP_H

#include "transform_base.h"

namespace ofec::free_peaks {
	class MapObjective : public TransformBase {
	private:
		std::vector<std::pair<Real, Real>> m_to, m_from;

	public:
		MapObjective(Problem *pro, const std::string &subspace_name, const ParameterMap &param);
		void transfer(std::vector<Real> &obj, const std::vector<Real>& var) override;
		void bindData() override;
	};
}

#endif // !OFEC_FREE_PEAKS_LINEAR_MAP_H
