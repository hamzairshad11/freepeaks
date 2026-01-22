#include "map_x_base.h"
#include "../../free_peaks.h"

namespace ofec::free_peaks {
	MapXBase::MapXBase(Problem *pro, const std::string& subspace_name, const ParameterMap& param) :
		TransformBase(pro, subspace_name, param), m_subspace_name(subspace_name){}

	void MapXBase::bindData() {


		m_subspace_var_ranges = CAST_FPs(m_pro)->subspaceTree().tree->getBox(m_subspace_name);
		auto& subpro = CAST_FPs(m_pro)->subspaceTree().name_box_subproblem.at(m_subspace_name).second;
		m_function_var_ranges = subpro->function()->varRanges();
		m_optimal_vars_inFunction = subpro->function()->optimalVarsFunction();
	
	}
}