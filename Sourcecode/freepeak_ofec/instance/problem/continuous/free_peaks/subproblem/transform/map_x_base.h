#ifndef OFEC_FREE_PEAKS_X_BASE_H
#define OFEC_FREE_PEAKS_X_BASE_H

#include "transform_base.h"
#include "../../../../../../core/problem/encoding.h"

namespace ofec::free_peaks {
	class MapXBase : public TransformBase {
	protected:
		const std::string m_subspace_name;
		std::vector<std::pair<Real, Real>> m_subspace_var_ranges;
		std::vector<std::pair<Real, Real>> m_function_var_ranges;
		std::vector<VariableVector<Real>> m_optimal_vars_inFunction;

	public:
		MapXBase(Problem *pro, const std::string& subspace_name, const ParameterMap& param);
		virtual void bindData() override;
	};
}

#endif // !OFEC_FREE_PEAKS_LINEAR_MAP_H
