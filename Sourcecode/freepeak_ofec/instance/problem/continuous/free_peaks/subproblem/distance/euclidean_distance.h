#ifndef OFEC_FREE_PEAKS_EUCLIDEAN_DISTANCE_H
#define OFEC_FREE_PEAKS_EUCLIDEAN_DISTANCE_H

#include "distance_base.h"

namespace ofec::free_peaks {
	class EuclideanDistance : public DistanceBase {
	public:
		EuclideanDistance(ofec::Problem *pro, const std::string &subspace_name, const ParameterMap &param);
		Real operator()(const std::vector<Real> &a, const std::vector<Real> &b) const override;
	};
}

#endif // !OFEC_FREE_PEAKS_EUCLIDEAN_DISTANCE_H
