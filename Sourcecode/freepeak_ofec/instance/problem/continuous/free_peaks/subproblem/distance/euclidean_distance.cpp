#include "euclidean_distance.h"

namespace ofec::free_peaks {


	void EuclideanDistance::initialize(Problem *pro, const std::string &subspace_name, const ParameterMap &param)  {
		DistanceBase::initialize(pro, subspace_name, param);
	}

	Real EuclideanDistance::operator()(const std::vector<Real> &a, const std::vector<Real> &b) const {
		if (a.size() != b.size())
			return ms_infinity;
		Real dis = 0;
		for (size_t i = 0; i < a.size(); ++i)
			dis += (a[i] - b[i]) * (a[i] - b[i]);
		return sqrt(dis);
	}
}