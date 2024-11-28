#pragma once

#include <array>

namespace fem::two_dim {

    struct MeshQuadLinear {
        std::array<size_t, 4> indOfPoints{}; /// indices of point in order [ BL BR TL RT ]
        size_t materialNum{}; /// Number of material of mesh
    };

}
