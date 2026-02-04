#include "transform_base.h"
#include "../../free_peaks.h"


namespace ofec::free_peaks {
	void TransformBase::initialize(Problem *pro, const std::string &subspace_name, const ParameterMap &param) {
		m_pro = pro;
		m_subspace_name = subspace_name;
		m_input_parameters.input(param);	
	}	

	void TransformBase::bindData() {
		m_sub_pro = CAST_FPs(m_pro)->subproblem(m_subspace_name);
	}
}
