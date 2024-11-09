#pragma once
#include <vector>

#include "two_dim_point.h"
#include "two_dim_quads_linear_mesh.h"
#include "two_dim_domain.h"

namespace fem::two_dim {

    struct GridQuadLinear {
        std::vector<Point> points; /// Points of Grid
        std::vector<MeshQuadLinear> meshes; /// Meshes of Grid
        std::vector<size_t> usedMaterials; /// Numbers of used materials for this grid

        void buildFrom(const Domain& domain);
    };

}
