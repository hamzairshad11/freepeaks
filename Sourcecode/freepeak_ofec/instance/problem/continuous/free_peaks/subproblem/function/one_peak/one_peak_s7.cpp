#include "one_peak_s7.h"
#include "../../../free_peaks.h"
#include "../../function/one_peak_function.h"

namespace ofec::free_peaks {
	OnePeakS7::OnePeakS7(Problem *pro, const std::string &subspace_name, const ParameterMap &param) :
		OnePeakBase(pro, subspace_name, param),
		m_radius(0) {}

	Real OnePeakS7::evaluate_(Real dummy, size_t var_size) {
		if (dummy <= m_radius)
			return m_height * cos(OFEC_PI * dummy / m_radius);
		else
			return -m_height - dummy + m_radius;
	}

	void OnePeakS7::bindData() {
		OnePeakBase::bindData();
		auto fun = dynamic_cast<OnePeakFunction*>(CAST_FPs(m_pro)->subproblem(m_subspace_name)->function());
		Real min_dis = fun->computeMinDis(m_center);
		if (m_param.has("r_ratio"))
			m_radius = m_param.get<Real>("r_ratio") * min_dis;
		else
			m_radius = 0.5 * min_dis;
	}
}