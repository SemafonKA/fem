#pragma once

#include <vector>
#include <array>
#include <functional>

#include "../../three_steppers/Headers/SparseMatrix.h";
#include "three_dim_cuboid_linear_grid.h"

namespace fem::three_dim {

    using LocalCuboidLinearMat = std::array<std::array<double, 8>, 8>;
    using LocalCuboidLinearVec = std::array<double, 8>;

    struct DomainFunctionsCuboidLinear {
        static constexpr double not_defined = std::numeric_limits<double>::infinity();

        std::function<double(double, double, double, size_t)> lambda; // function like double(x, y, z, material)
        std::function<double(double, double, double, size_t)> gamma;  // function like double(x, y, z, material)
        std::function<double(double, double, double, size_t)> func;   // function like double(x, y, z, material)
        std::function<double(double, double, double, size_t)> s1_func; // function like double(x, y, z, index). Must return DomainFunctionsCuboidLinear::not_defined if not defined for this coord line
        std::function<double(double, double, double, size_t)> s2_func; // function like double(x, y, z, index). Must return DomainFunctionsCuboidLinear::not_defined if not defined for this coord line
    };


    struct DomainFunctionsCuboidLinearDynamic {
        static constexpr double not_defined = std::numeric_limits<double>::infinity();

        std::function<double(double, double, double, size_t)> lambda; // function like double(x, y, z, material)
        std::function<double(double, double, double, size_t)> sigma;  // function like double(x, y, z, material)
        std::function<double(double, double, double, double, size_t)> func;    // function like double(x, y, z, t, material)
        std::function<double(double, double, double, double, size_t)> s1_func; // function like double(x, y, z, t, index). Must return DomainFunctionsCuboidLinear::not_defined if not defined for this coord line
        std::function<double(double, double, double, double, size_t)> s2_func; // function like double(x, y, z, t, index). Must return DomainFunctionsCuboidLinear::not_defined if not defined for this coord line
        std::function<double(double, double, double, double)> initials;        // function like double(x, y, z, t)
    };

    class SolverCuboidLinear {
    public:

        inline SolverCuboidLinear(const DomainFunctionsCuboidLinear& funcs, const GridCuboidLinear& grid) :
            _funcs(funcs), _grid(grid) {}

        [[nodiscard]]
        auto solveStatic() -> std::vector<double>;

        [[nodiscard]]
        auto value(Point p) -> double;

        inline auto getGrid() const -> const GridCuboidLinear& {
            return _grid;
        }

        inline auto getSolution() const -> const std::vector<double>& {
            return _solve;
        }

    private:
        SparseMatrix _global_mat;
        SparseMatrix _global_M;
        SparseMatrix _global_G;
        std::vector<double> _global_b;

        std::vector<double> _solve;

        DomainFunctionsCuboidLinear _funcs;
        GridCuboidLinear _grid;

        void generatePortrait();

        auto getLocalG(const MeshCuboidLinear& mesh) -> LocalCuboidLinearMat;
        auto getLocalM(const MeshCuboidLinear& mesh, bool addGamma = true) -> LocalCuboidLinearMat;
        auto getLocalB(const MeshCuboidLinear& mesh) -> LocalCuboidLinearVec;

        void includeS2();
        void includeS1();

        auto findFinite(const Point& p) const->size_t;
    };
}
