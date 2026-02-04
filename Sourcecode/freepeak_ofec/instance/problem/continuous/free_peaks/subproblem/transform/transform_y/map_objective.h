#ifndef OFEC_FREE_PEAKS_LINEAR_MAP_H
#define OFEC_FREE_PEAKS_LINEAR_MAP_H

#include "transform_y_base.h"

namespace ofec::free_peaks {
	class MapObjective : public Y_TransformBase {
		OFEC_CONCRETE_INSTANCE(MapObjective)
	private:
		std::vector<std::pair<Real, Real>> m_to, m_from;

	public:
		void addInputParameters() {}
		virtual void initialize(Problem *pro, const std::string &subspace_name, const ParameterMap &param)override;
		void transfer(std::vector<Real> &obj, const std::vector<Real>& var) override;
		virtual void bindData() override;
	};
}

#endif // !OFEC_FREE_PEAKS_LINEAR_MAP_H
