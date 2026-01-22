#include "register_transform.h"
#include "../../factory.h"
#include "map_objective.h"
//#include "map_Weierstrass_rugeness.h"
#include "map_x_irregularities.h"
#include "map_x_assymetrix.h"
#include "map_x_ill_conditioning.h"
namespace ofec::free_peaks {
	void registerTransform() {
		REGISTER_FP(TransformBase, MapObjective, "map_objective");
		//REGISTER_FP(TransformBase, WeierstrassRugenessMap, "WeierstrassRugenessMap");
		REGISTER_FP(TransformBase, MapXIrregularity, "MapXIrregularity");
		REGISTER_FP(TransformBase, MapXAssymetrix, "MapXAssymetrix");
		REGISTER_FP(TransformBase, MapXIllConditioning, "MapXIllConditioning");
	}
}