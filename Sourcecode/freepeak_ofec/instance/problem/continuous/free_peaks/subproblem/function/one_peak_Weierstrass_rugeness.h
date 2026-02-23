#ifndef OFEC_FREE_PEAKS_WEIERSTRASS_RUG_H
#define OFEC_FREE_PEAKS_WEIERSTRASS_RUG_H

#include "transform_x_base.h"
#include "../../../../../../../instance/problem/continuous/single_objective/global/bbob/bbob.h"


// to add

namespace ofec::free_peaks {
	class WeierstrassRugenessMap : public OnePeakBase, public BBOB_base {
		OFEC_CONCRETE_INSTANCE(WeierstrassRugenessMap)
	private:

		Real m_alpha = 0.5;        // a = 0.5
		int m_b = 3.0;            // b = 3
		int m_kmax = 20;          // CEC2005 标准是 20 项

	public:
		void addInputParameters();
		virtual void initialize(Problem *pro, const std::string& subspace_name, const ParameterMap& param)override;
		void transfer(std::vector<Real>& obj, const std::vector<Real>& var) override;
		virtual Real evaluate_(Real dummy, size_t var_size)override;
		
	 	virtual void bindData()override;
	};
}

#endif // !OFEC_FREE_PEAKS_LINEAR_MAP_H
