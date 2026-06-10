#include "map_objective.h"
#include "../../../free_peaks.h"

namespace ofec::free_peaks {

	void MapObjective::initialize(Problem* pro, const std::string& subspace_name, const ParameterMap& param) {
		TransformBase::initialize(pro, subspace_name, param);
		//bindData();
	}

	
	void MapObjective::transfer(std::vector<Real> &obj, const std::vector<Real>& var) {
		for (size_t j = 0; j < m_pro->numberObjectives(); ++j) {
			if (m_from[j].second != m_from[j].first) {
				obj[j] = m_to[j].first + (m_to[j].second - m_to[j].first)
					* (obj[j] - m_from[j].first) / (m_from[j].second - m_from[j].first);
				if (obj[j] < m_to[j].first) obj[j] = m_to[j].first;
				else if (obj[j] > m_to[j].second) obj[j] = m_to[j].second;
			}
		}
	}

	void MapObjective::bindData() {
		Y_TransformBase::bindData();
		
		auto* subpro = CAST_FPs(m_pro)->subproblem(m_subspace_name);
		m_to = m_from = subpro->function()->objRanges();
		for (size_t j = 0; j < m_to.size(); ++j) {
			 m_to[j].first = 0;
			 m_to[j].second = 1;
		}
	}
}