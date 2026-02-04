#ifndef OFEC_FREE_PEAKS_FLAT_BORDER_H
#define OFEC_FREE_PEAKS_FLAT_BORDER_H

#include "distance_base.h"
#include "../function/function_base.h"

namespace ofec::free_peaks {
	class FlatBorder : public DistanceBase {
		OFEC_CONCRETE_INSTANCE(FlatBorder)

	private:
		FunctionBase *m_fun;

	public:

		void addInputParameters() {}
		virtual void initialize(Problem *pro, const std::string &subspace_name, const ParameterMap &param) override;
		Real operator()(const std::vector<Real> &a, const std::vector<Real> &b) const override;
		void bindData() override;
	};
}

#endif // !OFEC_FREE_PEAKS_FLAT_BORDER_H
