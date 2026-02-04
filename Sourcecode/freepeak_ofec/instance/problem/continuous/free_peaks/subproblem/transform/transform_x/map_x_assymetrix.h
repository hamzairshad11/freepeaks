#ifndef OFEC_FREE_PEAKS_X_ASSYMETRIX_H
#define OFEC_FREE_PEAKS_X_ASSYMETRIX_H

#include "transform_x_base.h"
#include "map_x_base.h"

namespace ofec::free_peaks {
	class MapXAssymetrix : public MapXBase {
		OFEC_CONCRETE_INSTANCE(MapXAssymetrix)
	private:
	public:

		void addInputParameters() {}	
		virtual void initialize(Problem *pro, const std::string& subspace_name, const ParameterMap& param)override;
		void transfer(std::vector<Real>& x, const std::vector<Real>& var) override;

	};
}

#endif // !OFEC_FREE_PEAKS_LINEAR_MAP_H
