#include "map_x_irregularities.h"
#include "../../../free_peaks.h"
#include "../../../../../../../utility/functional.h"

namespace ofec::free_peaks {

	MapXIrregularity* MapXIrregularity::create() {
		return new MapXIrregularity();
	}

	void MapXIrregularity::initialize(Problem* pro, const std::string& subspace_name, const ParameterMap& param) {
		TransformBase::initialize(pro, subspace_name, param);
	}

	void MapXIrregularity::transfer(std::vector<Real>& x, const std::vector<Real>& var) {
		double c1, c2, a;
		for (int i = 0; i < x.size(); ++i) {
			if (x[i] > 0) {
				c1 = 10;	c2 = 7.9;
			}
			else {
				c1 = 5.5;	c2 = 3.1;
			}
			if (x[i] != 0) {
				a = log(fabs(x[i]));
			}
			else a = 0;
			x[i] = sign(x[i]) * exp(a + 0.049 * (sin(c1 * x[i] + sin(c2 * x[i]))));
		}
	}

}