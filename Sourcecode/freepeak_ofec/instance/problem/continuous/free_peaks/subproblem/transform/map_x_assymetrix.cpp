#include "map_x_assymetrix.h"
#include "../../free_peaks.h"

namespace ofec::free_peaks {
	MapXAssymetrix::MapXAssymetrix(Problem *pro, const std::string& subspace_name, const ParameterMap& param) :
		MapXBase(pro, subspace_name, param) {}

	void MapXAssymetrix::transfer(std::vector<Real>& x_, const std::vector<Real>& var) {


		if (x_.size()==1) return;

		double belta = 0.1;
		//for (int idx(0); idx < x_.size(); ++idx) {
		//	x_[idx] = ofec::mapReal(x_[idx], m_subspace_var_ranges[idx].first, m_subspace_var_ranges[idx].second,
		//		m_function_var_ranges[idx].first, m_function_var_ranges[idx].second
		//	);
		//}


		//for (int idx(0); idx < x_.size(); ++idx) {
		//	x_[idx] -= m_optimal_vars_inFunction[idx];
		//}

		for (int i = 0; i < x_.size(); ++i) {
			if (x_[i] > 0) {
				x_[i] = pow(x_[i], 1. + belta * (i / (x_.size() - 1)) * log(x_[i] + 1));
			}
		}

		//for (int idx(0); idx < x_.size(); ++idx) {
		//	x_[idx] += m_optimal_vars_inFunction.front()[idx];
		//}

	
	}

	void MapXAssymetrix::bindData() {
		MapXBase::bindData();
		
	}
}