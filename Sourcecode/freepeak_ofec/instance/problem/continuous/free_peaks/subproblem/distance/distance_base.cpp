#include "distance_base.h"
#include <limits>
#include "../../../../../../core/problem/problem.h"

namespace ofec::free_peaks {
	const Real DistanceBase::ms_infinity = std::numeric_limits<Real>::max();

	DistanceBase::DistanceBase(ofec::Problem *pro, const std::string &subspace_name, const ParameterMap &param) : 
		m_pro(pro), m_subspace_name(subspace_name), m_param(param) {}
}
