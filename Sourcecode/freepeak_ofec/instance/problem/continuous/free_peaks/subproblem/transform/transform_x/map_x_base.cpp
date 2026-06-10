#include "map_x_base.h"
#include "../../../free_peaks.h"

namespace ofec::free_peaks {
	void MapXBase::initialize(Problem* pro, const std::string& subspace_name, const ParameterMap& param) {
		TransformBase::initialize(pro, subspace_name, param);
	}

	void MapXBase::bindData() {

		X_TransformBase::bindData();

		m_subspace_var_ranges = CAST_FPs(m_pro)->subspaceBox(m_subspace_name);
		auto* subpro = CAST_FPs(m_pro)->subproblem(m_subspace_name);
		m_function_var_ranges = subpro->function()->varRanges();
		m_optimal_vars_inFunction = subpro->function()->optimalVarsFunction();
	
	}
}