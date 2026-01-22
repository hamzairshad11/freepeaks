#include "init_pop_bounded.h"

namespace ofec {
	void InitPopBounded::initializeVariables(VariableBase &x_, Random *rnd) const {
		auto &x = dynamic_cast<VariableVector<Real>&>(x_);
		for (int i = 0; i < m_number_variables; ++i) {
			if (m_domain_init_pop[i].limited)
				x[i] = rnd->uniform.nextNonStd(m_domain_init_pop[i].limit.first, m_domain_init_pop[i].limit.second);
			else {
				if (m_domain[i].limited)     // If m_initial_domain is not given, then use problem boundary as initialization range
					x[i] = rnd->uniform.nextNonStd(m_domain[i].limit.first, m_domain[i].limit.second);
				else                         // Else if the problem function has no boundary 
					x[i] = rnd->uniform.nextNonStd(-std::numeric_limits<Real>::max(), std::numeric_limits<Real>::max());
			}
		}
	}

	void InitPopBounded::resizeVariable(size_t num_vars) {
		Continuous::resizeVariable(num_vars);
		m_domain_init_pop.resize(num_vars);
	}

	void InitPopBounded::setDomainInitPop(Real l, Real u) {
		for (size_t i = 0; i < m_number_variables; ++i)
			m_domain_init_pop.setRange(l, u, i);
	}

	void InitPopBounded::setDomainInitPop(const Domain<Real> &domain) {
		m_domain_init_pop = domain;
	}

	void InitPopBounded::clearDomainInitPop() {
		for (size_t i = 0; i < m_number_variables; ++i)
			m_domain_init_pop[i].limited = false;
	}

	void InitPopBounded::setDomainInitPop(const std::vector<std::pair<Real, Real>> &r) {
		size_t count = 0;
		for (auto &i : r) {
			m_domain_init_pop.setRange(i.first, i.second, count++);
		}
	}

	void InitPopBounded::initialize_(Environment *env) {
		Continuous::initialize_(env);
		m_domain_init_pop.resize(m_number_variables);
	}
}