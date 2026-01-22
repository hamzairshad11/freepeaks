#include"constraint_base.h"
#include "../../free_peaks.h"

namespace ofec::free_peaks {
	ConstraintBase::ConstraintBase(Problem *pro, const std::string &subspace_name, const ParameterMap &param) :
		m_pro(pro), m_subspace_name(subspace_name), m_param(param), m_dis(nullptr) {}

	void ConstraintBase::bindData() {
		m_dis = CAST_FPs(m_pro)->subproblem(m_subspace_name)->distance();
	}
}
