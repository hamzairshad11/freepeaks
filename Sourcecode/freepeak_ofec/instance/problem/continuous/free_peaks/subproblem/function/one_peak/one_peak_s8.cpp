#include "one_peak_s8.h"
#include "../../../free_peaks.h"
#include "../../function/one_peak_function.h"

namespace ofec::free_peaks {

	void OnePeakS8::addInputParameters() {
		m_input_parameters.add("m", new RangedInt(m_m, 0, 100, 3));
		m_input_parameters.add("eta", new RangedReal(m_eta, 0, 1e6, 5.5));
		m_input_parameters.add("r_ratio", new RangedReal(m_r_ratio, 0, 1, 0.5));
	}

	void OnePeakS8::initialize(Problem* pro, const std::string& subspace_name, const ParameterMap& param) {
		OnePeakBase::initialize(pro, subspace_name, param);
	}

	Real OnePeakS8::evaluate_(Real dummy, size_t var_size) {
		// Guard: if radius is zero (from degenerate transforms), return peak height.
		if (m_radius <= Real(1e-12))
			return m_height;
		if (dummy <= m_radius) {
			int n = static_cast<int>(floor(m_m * dummy / m_radius));
			return m_height * (cos(m_m * OFEC_PI * (dummy - n * m_radius / m_m) / m_radius)
				- m_eta * n) / (sqrt(dummy + 1));
		}
		else {
			return m_height * (-1 - m_eta * (m_m - 1)) / (sqrt(m_radius + 1)) - (dummy - m_radius);
		}
	}

	void OnePeakS8::bindData() {
		OnePeakBase::bindData();
		auto fun = dynamic_cast<OnePeakFunction*>(CAST_FPs(m_pro)->subproblem(m_subspace_name)->function());
		Real min_dis = fun->computeMinDis(m_center);
		m_radius = m_r_ratio * min_dis;

	}
}