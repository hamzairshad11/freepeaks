#include "one_peak_s7.h"
#include "../../../free_peaks.h"
#include "../../function/one_peak_function.h"

namespace ofec::free_peaks {

	void OnePeakS7::addInputParameters() {
		m_input_parameters.add("r_ratio", new RangedReal(m_r_ratio, 0, 1, 0.5));
	}



	void OnePeakS7::initialize(Problem* pro, const std::string& subspace_name, const ParameterMap& param) {
		OnePeakBase::initialize(pro, subspace_name, param);
		//	bindData();
	}

	Real OnePeakS7::evaluate_(Real dummy, size_t var_size) {
		// Guard: if radius is zero or negative (from degenerate transforms),
		// return the peak height to avoid cos(pi/0) = NaN.
		if (m_radius <= Real(1e-12))
			return m_height;
		if (dummy <= m_radius)
			return m_height * cos(OFEC_PI * dummy / m_radius);
		else
			return -m_height - dummy + m_radius;
	}

	void OnePeakS7::bindData() {
		OnePeakBase::bindData();
		auto fun = dynamic_cast<OnePeakFunction*>(CAST_FPs(m_pro)->subproblem(m_subspace_name)->function());
		Real min_dis = fun->computeMinDis(m_center);
		m_radius = m_r_ratio * min_dis;
	}
}