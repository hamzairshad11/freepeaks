#ifndef OFEC_MAP_X_ROTATION_H
#define OFEC_MAP_X_ROTATION_H

#include "map_x_base.h"

namespace ofec::free_peaks {

    class MapXRotation : public MapXBase {
        OFEC_CONCRETE_INSTANCE(MapXRotation)
    public:
        void initialize(Problem* pro, const std::string& subspace_name, const ParameterMap& param) override;
        void transfer(std::vector<Real>& x_, const std::vector<Real>& var) override;

    private:
        std::vector<std::vector<Real>> m_rotation_matrix;   // orthogonal matrix R (dim x dim)

        void generateRandomRotationMatrix(size_t dim, unsigned seed);
        void gramSchmidt(std::vector<std::vector<Real>>& matrix);
        Real computeDeterminant(const std::vector<std::vector<Real>>& matrix);
    };

} // namespace ofec::free_peaks

#endif // OFEC_MAP_X_ROTATION_H