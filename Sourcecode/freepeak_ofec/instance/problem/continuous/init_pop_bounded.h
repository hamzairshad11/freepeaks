#ifndef OFEC_INIT_POP_BOUNDED_H
#define OFEC_INIT_POP_BOUNDED_H

#include "../../../core/problem/continuous/continuous.h"

namespace ofec {
	class InitPopBounded : virtual public Continuous {
	public:
		void initializeVariables(VariableBase &x, Random *rnd) const override;
		void resizeVariable(size_t num_vars) override;
		const std::pair<Real, Real>& rangeInitPop(size_t i) const { return m_domain_init_pop.range(i).limit; }
		const Domain<Real>& domainInitPop() const { return m_domain_init_pop; }
		void setDomainInitPop(Real l, Real u);
		void setDomainInitPop(const std::vector<std::pair<Real, Real>> &r);
		void setDomainInitPop(const Domain<Real> &domain);
		void clearDomainInitPop();

	protected:
		void initialize_(Environment *env) override;

		Domain<Real> m_domain_init_pop; // domain for intializing population
	};
}

#endif // !OFEC_INIT_POP_BOUNDED_H
