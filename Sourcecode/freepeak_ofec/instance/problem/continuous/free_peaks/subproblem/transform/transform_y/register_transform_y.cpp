#include "register_transform_y.h"
#include "../../../factory.h"
#include "map_objective.h"

namespace ofec::free_peaks {
	void registerTransformY() {
		REGISTER_FP(Y_TransformBase, MapObjective, "map_objective");
	}
}