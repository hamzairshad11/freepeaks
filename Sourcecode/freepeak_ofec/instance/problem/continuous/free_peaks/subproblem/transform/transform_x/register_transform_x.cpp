#include "register_transform_x.h"
#include "../../../factory.h"
//#include "map_objective.h"
//#include "map_Weierstrass_rugeness.h"
#include "map_x_irregularities.h"
#include "map_x_assymetrix.h"
#include "map_x_ill_conditioning.h"
#include "map_x_partybias.h"
#include "map_x_deceptive.h"
#include "map_x_linkage.h"

namespace ofec::free_peaks {
	void registerTransformX() {
		REGISTER_FP(X_TransformBase, MapXIrregularity, "MapXIrregularity");
		REGISTER_FP(X_TransformBase, MapXAssymetrix, "MapXAssymetrix");
		REGISTER_FP(X_TransformBase, MapXIllConditioning, "MapXIllConditioning");
		REGISTER_FP(X_TransformBase, MapXPartyBias, "MapXPartyBias");
		REGISTER_FP(X_TransformBase, MapXDeceptive, "MapXDeceptive");
		REGISTER_FP(X_TransformBase, MapXLinkage, "MapXLinkage");
	}
}