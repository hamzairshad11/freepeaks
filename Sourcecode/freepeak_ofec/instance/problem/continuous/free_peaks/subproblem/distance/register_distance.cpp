#include "register_distance.h"
#include "../../factory.h"
#include "euclidean_distance.h"
#include "flat_border.h"

namespace ofec::free_peaks {
	void registerDistance() {
		REGISTER_FP(DistanceBase, EuclideanDistance, "Euclidean");
		REGISTER_FP(DistanceBase, FlatBorder, "flat_border");
	}
}