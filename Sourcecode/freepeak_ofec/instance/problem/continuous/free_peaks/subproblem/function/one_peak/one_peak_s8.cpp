#include "one_peak_s8.h"
#include "../../../free_peaks.h"
#include "../../function/one_peak_function.h"

namespace ofec::free_peaks {
	OnePeakS8::OnePeakS8(Problem *pro, const std::string &subspace_name, const ParameterMap &param) :
		OnePeakBase(pro, subspace_name, param) ,
		m_radius(0),
		m_m(param.get<int>("m")),
		m_eta(param.get<Real>("eta")) {}

	Real OnePeakS8::evaluate_(Real dummy, size_t var_size) {
		if (dummy <= m_radius) {
			int n = floor(m_m * dummy / m_radius);
			return m_height *(cos(m_m * OFEC_PI * (dummy - n * m_radius / m_m) / m_radius) - m_eta * n) / (sqrt(dummy + 1));
		}
		else {
			return m_height * (-1 - m_eta * (m_m - 1)) / (sqrt(m_radius + 1)) - (dummy - m_radius);
		}
	}

	void OnePeakS8::bindData() {
		OnePeakBase::bindData();
		auto fun = dynamic_cast<OnePeakFunction*>(CAST_FPs(m_pro)->subproblem(m_subspace_name)->function());
		Real min_dis = fun->computeMinDis(m_center);
		if (m_param.has("r_ratio"))
			m_radius = m_param.get<Real>("r_ratio") * min_dis;
		else
			m_radius = 0.5 * min_dis;
	}
}