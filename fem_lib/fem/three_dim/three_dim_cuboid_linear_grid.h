#pragma once
#include <vector>

#include "three_dim_point.h"
#include "three_dim_cuboid_linear_mesh.h"
#include "three_dim_domain.h"

namespace fem::three_dim {

    struct GridCuboidLinear {
        size_t Kx{}; /// Count of points by X axis
        size_t Ky{}; /// Count of points by Y axis
        size_t Kz{}; /// Count of points by Z axis

        std::vector<Point> points{};          /// Points of Grid
        std::vector<MeshCuboidLinear> meshes{}; /// Meshes of Grid
        std::vector<size_t> usedMaterials{};  /// Numbers of used materials for this grid

        /**
         * @brief Build grid from provided domain
         * @param domain - the domain from which the grid will be built
         */
        void buildFrom(const Domain& domain);

        /**
        * @brief Format grid points to string in format like `x[0] y[0] z[0]   x[1] y[1] z[1]   ... \\n x[Kx] y[Kx] z[Kx]   ...`
        * @return formatted string
        */
        [[nodiscard]]
        auto dumpPoints() const->std::string;

        /**
         * @brief Format grid meshes to string in format like `m[i]   p[0] p[1] p[2] p[3] p[4] p[5] p[6] p[7] \n ...`
         * @return formatted string
         */
        [[nodiscard]]
        auto dumpMeshes() const->std::string;

        /**
         * @brief Dump grid to string if computer-readable format (see `dumpPoints` and `dumpMeshes` to know what format is)
         * @return formatted string
         */
        [[nodiscard]]
        auto dump() const->std::string;
    };

}
