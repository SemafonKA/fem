#pragma once

#include <array>
#include <vector>
#include <format>
#include <sstream>

namespace fem::three_dim {

    struct EdgeCondition {
        enum class Type { First, Second };
        Type type{};
        size_t functionNumber{ 0 };
        std::array<size_t, 4> pointsIndices{}; /// indices of two connected points 

        inline auto toString() const -> std::string {
            return std::format(
                "[ type: {}, funNum: {}, points: [{:5}, {:5}, {:5}, {:5}] ]",
                (type == Type::First) ? "First" : "Second",
                functionNumber,
                pointsIndices[0],
                pointsIndices[1],
                pointsIndices[2],
                pointsIndices[3]
            );
        }
    };

    struct MeshQuadLinear {
        std::array<size_t, 12> indOfEdges{};  /// indices of edges in order [  ]
        std::array<size_t, 8> indOfPoints{}; /// indices of point in order [ BFL BFR BRL BRR TFL TFT TRL TRR ] ([B]ottom, [T]op, [L]eft, [R]ight, [F]ront, [R]ear)
        size_t materialNum{}; /// Number of material of mesh
        std::vector<EdgeCondition> s1Conditions{};  /// s1 Edge conditions
        std::vector<EdgeCondition> s2Conditions{}; /// s2 Edge conditions

        inline auto toString() const -> std::string {
            std::stringstream ss;
            ss << "[ material: " << materialNum << ", points: [ " << indOfPoints[0];
            for (size_t i = 1; i < indOfPoints.size(); i++) {
                ss << ", " << indOfPoints[i];
            }
            ss << " ], edges: [ " << indOfEdges[0];
            for (size_t i = 1; i < indOfEdges.size(); i++) {
                ss << ", " << indOfEdges[i];
            }
            ss << " ], edgeConditions: {";
            for (const auto& edgeCond : s1Conditions) {
                ss << "\n  " << edgeCond.toString();
            }
            for (const auto& edgeCond : s2Conditions) {
                ss << "\n  " << edgeCond.toString();
            }
            ss << "\n}]";
            return ss.str();
        }
    };

}
