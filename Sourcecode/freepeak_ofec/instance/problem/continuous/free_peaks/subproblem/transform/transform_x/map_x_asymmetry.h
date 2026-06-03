#ifndef OFEC_MAP_X_ASYMMETRY_H
#define OFEC_MAP_X_ASYMMETRY_H

#include "map_x_base.h"

namespace ofec::free_peaks {

    class MapXAsymmetry : public MapXBase {
		OFEC_CONCRETE_INSTANCE(MapXAsymmetry)
    public:
        void initialize(Problem* pro, const std::string& subspace_name, const ParameterMap& param) override;
        void transfer(std::vector<Real>& x_, const std::vector<Real>& var) override;

    private:
        Real m_beta = Real(0.2);
    };

} // namespace ofec::free_peaks

#endif // OFEC_MAP_X_ASYMMETRY_H