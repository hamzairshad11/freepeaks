#ifndef OFEC_FREE_PEAKS_X_IRREGULARITIES_H
#define OFEC_FREE_PEAKS_X_IRREGULARITIES_H

#include "transform_base.h"

namespace ofec::free_peaks {
	class MapXIrregularity : public TransformBase {
	private:

	public:
		MapXIrregularity(Problem *pro, const std::string& subspace_name, const ParameterMap& param);
		void transfer(std::vector<Real>& x, const std::vector<Real>& var) override;
		void bindData() override;
	};
}

#endif // !OFEC_FREE_PEAKS_LINEAR_MAP_H
