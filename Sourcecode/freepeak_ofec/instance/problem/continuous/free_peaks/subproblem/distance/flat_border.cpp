#include "flat_border.h"
#include "../../free_peaks.h"
#include "../transform/map_objective.h"

namespace ofec::free_peaks {
	FlatBorder::FlatBorder(Problem *pro, const std::string &subspace_name, const ParameterMap &param) :
		DistanceBase(pro, subspace_name, param), m_fun(nullptr) {}

	Real FlatBorder::operator()(const std::vector<Real> &a, const std::vector<Real> &b) const {
		if (a.size() != b.size())
			return ms_infinity;
		auto &ranges = m_fun->varRanges();
		Real dis = 1.0, val;
		for (size_t i = 0; i < a.size(); ++i) {
			val = a[i] - b[i];
			if (val >= 0)
				val = 1 - val / (ranges[i].second - b[i]);
			else
				val = 1 + val / (b[i] - ranges[i].first);
			dis *= val;
		}
		dis = (1 - dis) * m_fun->domainSize();
		return dis;
	}
	
	void FlatBorder::bindData() {
		m_fun = CAST_FPs(m_pro)->subproblem(m_subspace_name)->function();
	}
}