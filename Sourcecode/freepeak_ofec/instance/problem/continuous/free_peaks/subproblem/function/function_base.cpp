#include"function_base.h"
#include "../../free_peaks.h"

namespace ofec::free_peaks {
	FunctionBase::FunctionBase(Problem *pro, const std::string &subspace_name, const ParameterMap &param) : 
		m_pro(pro), m_subspace_name(subspace_name), m_param(param), m_domain_size(0),
		m_optimal_vars_mapped(true) {}

	void FunctionBase::evaluate(const std::vector<Real> &var, std::vector<double> &obj) {
		std::vector<Real> var_(var.size());
		auto& from = CAST_FPs(m_pro)->subspaceTree().tree->getBox(m_subspace_name);
		auto& to = m_var_ranges;
		for (size_t i = 0; i < var.size(); ++i) {
			var_[i] = (var[i] - from[i].first) / (from[i].second - from[i].first) *
				(to[i].second - to[i].first) + to[i].first;
		}
	//	transferX(var_, var);
		evaluate_(var_, obj);
	}

	void FunctionBase::transferX(std::vector<Real>& vecX, const std::vector<Real>& var) const {
		
		auto subpro = CAST_FPs(m_pro)->subproblem(m_subspace_name);
		return subpro->transferX(vecX, var);
	}
	
	void FunctionBase::mapOptimalVars() {
		if (!m_optimal_vars_mapped) {
			auto &from = CAST_FPs(m_pro)->subspaceTree().tree->getBox(m_subspace_name);
			auto &to = m_var_ranges;
			for (auto &var : m_optimal_vars) {
				for (size_t i = 0; i < var.size(); ++i) {
					var[i] = (var[i] - to[i].first) / (to[i].second - to[i].first) *
						(from[i].second - from[i].first) + from[i].first;
				}
			}
			m_optimal_vars_mapped = true;
		}
	}
}
