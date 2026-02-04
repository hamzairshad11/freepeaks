#ifndef OFEC_FREE_PEAKS_X_IRREGULARITIES_H
#define OFEC_FREE_PEAKS_X_IRREGULARITIES_H

#include "transform_x_base.h"
#include "../../../../../../../core/instance.h"

namespace ofec::free_peaks {
	class MapXIrregularity : public X_TransformBase {
		OFEC_CONCRETE_INSTANCE(MapXIrregularity)
	private:

	public:

		void addInputParameters() {}	

		virtual void initialize(Problem *pro, const std::string& subspace_name, const ParameterMap& param)override;
		void transfer(std::vector<Real>& x, const std::vector<Real>& var) override;

	};
}

#endif // !OFEC_FREE_PEAKS_LINEAR_MAP_H
