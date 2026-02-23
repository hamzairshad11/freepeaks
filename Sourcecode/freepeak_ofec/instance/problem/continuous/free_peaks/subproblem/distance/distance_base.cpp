#include "distance_base.h"
#include <limits>
#include "../../../../../../core/problem/problem.h"

namespace ofec::free_peaks {
	const Real DistanceBase::ms_infinity = std::numeric_limits<Real>::max();

	void DistanceBase::addInputParameters() {
		m_input_parameters.add("distance", new InputParameterValueTypeString(m_register_name));
	}
	void DistanceBase::initialize(ofec::Problem* pro, const std::string& subspace_name, const ParameterMap& param) {
		m_pro = (pro);
		m_subspace_name = (subspace_name);
		m_input_parameters.input(param);
		//	recordInputParameters();
	}
}