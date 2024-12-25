#pragma once

#include <array>
#include <vector>
#include <format>
#include <sstream>

namespace fem::two_dim {

    struct EdgeCondition {
        enum class Type { First, Second };
        Type type{};
        size_t functionNumber {0};
        std::array<size_t, 2> pointsIndices{}; /// indices of two connected points 

        inline auto toString() const -> std::string {
            return std::format(
                "[ type: {}, funNum: {}, points: [{:5}, {:5}] ]",
                (type == Type::First) ? "First" : "Second",
                functionNumber,
                pointsIndices[0],
                pointsIndices[1]
            );
        }
    };

    struct MeshQuadLinear {
        std::array<size_t, 4> indOfEdges{};  /// indices of edges in order [B L R T]
        std::array<size_t, 4> indOfPoints{}; /// indices of point in order [ BL BR TL RT ]
        size_t materialNum{}; /// Number of material of mesh
        std::vector<EdgeCondition> s1Conditions{}; /// s1 Edge conditions
        std::vector<EdgeCondition> s2Conditions{}; /// s1 Edge conditions

        inline auto toString() const -> std::string {
            std::stringstream ss;
            ss << "[ material: "
                << materialNum
                << ", points: ["
                << indOfPoints[0] << ", " << indOfPoints[1] << ", "
                << indOfPoints[2] << ", " << indOfPoints[3]
                << "], edges: ["
                << indOfEdges[0] << ", " << indOfEdges[1] << ", "
                << indOfEdges[2] << ", " << indOfEdges[3]
                << "], edgeConditions: {";
            for (const auto& edgeCond : s1Conditions) {
                ss << "\n  " << edgeCond.toString();
            }
            for (const auto& edgeCond : s2Conditions) {
                ss << "\n  " << edgeCond.toString();
            }
            ss << "\n}]";
            return ss.str();
            //return std::format(
            //    "[ material: {}, points: [{:5}, {:5}, {:5}, {:5}], edges: [{:5}, {:5}, {:5}, {:5}] ]",
            //    materialNum,
            //    indOfPoints[0],
            //    indOfPoints[1],
            //    indOfPoints[2],
            //    indOfPoints[3],
            //    indOfEdges[0],
            //    indOfEdges[1],
            //    indOfEdges[2],
            //    indOfEdges[3]
            //);
        }
    };

}
