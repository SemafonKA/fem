#pragma once
#include <vector>

#include "two_dim_point.h"
#include "two_dim_quads_linear_mesh.h"
#include "two_dim_domain.h"

namespace fem::two_dim {

    struct GridQuadLinear {
        size_t Kx{}; /// Count of points by X axis
        size_t Ky{}; /// Count of points by Y axis

        std::vector<Point> points{}; /// Points of Grid
        std::vector<MeshQuadLinear> meshes{}; /// Meshes of Grid
        std::vector<size_t> usedMaterials{}; /// Numbers of used materials for this grid

        void buildFrom(const Domain& domain);
    };

}
