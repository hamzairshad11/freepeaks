#include "flat_border.h"
#include "../../free_peaks.h"
#include <cmath>
#include "../transform/transform_y/map_objective.h"

namespace ofec::free_peaks {
	void FlatBorder::initialize(Problem* pro, const std::string& subspace_name, const ParameterMap& param) {
		DistanceBase::initialize(pro, subspace_name, param);
	}



	Real FlatBorder::operator()(const std::vector<Real>& a, const std::vector<Real>& b) const {
		if (a.size() != b.size())
			return ms_infinity;
		auto& ranges = m_fun->varRanges();
		Real dis = 1.0, val;
		for (size_t i = 0; i < a.size(); ++i) {
			val = a[i] - b[i];
			if (val >= 0) {
				const Real denom = ranges[i].second - b[i];
				val = std::abs(denom) > 1e-12 ? 1 - val / denom : (std::abs(val) <= 1e-12 ? 1 : 0);
			}
			else {
				const Real denom = b[i] - ranges[i].first;
				val = std::abs(denom) > 1e-12 ? 1 + val / denom : (std::abs(val) <= 1e-12 ? 1 : 0);
			}
			if (val < 0) val = 0;
			else if (val > 1) val = 1;
			dis *= val;
		}
		dis = (1 - dis) * m_fun->domainSize();
		return dis;
	}

	void FlatBorder::bindData() {
		DistanceBase::bindData();
		m_fun = CAST_FPs(m_pro)->subproblem(m_subspace_name)->function();
	}
}