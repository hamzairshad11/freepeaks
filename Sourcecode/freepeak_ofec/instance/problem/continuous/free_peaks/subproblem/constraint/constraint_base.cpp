#include"constraint_base.h"
#include "../../free_peaks.h"

namespace ofec::free_peaks {
	void ConstraintBase::initialize(Problem *pro, const std::string &subspace_name, const ParameterMap &param) {
		m_pro = pro;
		m_subspace_name = subspace_name;
		m_input_parameters.input(param);
	//	recordInputParameters();
	
	}

	void ConstraintBase::bindData() {
		m_dis = CAST_FPs(m_pro)->subproblem(m_subspace_name)->distance();
	}


}
