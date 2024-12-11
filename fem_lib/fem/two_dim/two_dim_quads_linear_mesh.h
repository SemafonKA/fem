#pragma once

#include <array>
#include <format>

namespace fem::two_dim {

    struct MeshQuadLinear {
        std::array<size_t, 4> indOfEdges{};  /// indices of edges in order [B L R T]
        std::array<size_t, 4> indOfPoints{}; /// indices of point in order [ BL BR TL RT ]
        size_t materialNum{}; /// Number of material of mesh

        inline auto toString() const -> std::string {
            return std::format(
                "[ material: {}, points: [{}, {}, {}, {}], edges: [{}, {}, {}, {}] ]",
                materialNum,
                indOfPoints[0],
                indOfPoints[1],
                indOfPoints[2],
                indOfPoints[3],
                indOfEdges[0],
                indOfEdges[1],
                indOfEdges[2],
                indOfEdges[3]
            );
        }
    };

}
