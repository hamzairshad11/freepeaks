#include "transform_base.h"

namespace ofec::free_peaks {
	TransformBase::TransformBase(Problem *pro, const std::string &subspace_name, const ParameterMap &param) :
		m_pro(pro), m_subspace_name(subspace_name), m_param(param) {}
}
