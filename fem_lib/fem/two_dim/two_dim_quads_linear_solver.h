#pragma once

#include <vector>
#include <array>
#include <functional>

#include "../../three_steppers/Headers/SparseMatrix.h";
#include "two_dim_quads_linear_grid.h"

namespace fem::two_dim {

    using LocalQuadsLinearMat = std::array<std::array<double, 4>, 4>;
    using LocalQuadsLinearVec = std::array<double, 4>;

    struct DomainFunctionsQuadsLinear {
        static constexpr double not_defined = std::numeric_limits<double>::infinity();

        std::function<double(double, double, size_t)> lambda; // function like double(x, y, material)
        std::function<double(double, double, size_t)> gamma;  // function like double(x, y, material)
        std::function<double(double, double, size_t)> func;   // function like double(x, y, material)
        std::function<double(double, double, size_t)> s1_func; // function like double(x, y, index). Must return DomainFunctionsQuadsLinear::not_defined if not defined for this coord line
        std::function<double(double, double, size_t)> s2_func; // function like double(x, y, index). Must return DomainFunctionsQuadsLinear::not_defined if not defined for this coord line
    };

    class SolverQuadsLinear {
    public:

        inline SolverQuadsLinear(const DomainFunctionsQuadsLinear& funcs, const GridQuadLinear& grid) :
            _funcs(funcs), _grid(grid) {}

        [[nodiscard]]
        auto solveStatic() -> std::vector<double>;

        [[nodiscard]]
        auto value(Point p) -> double;

    private:
        SparseMatrix _global_mat;
        SparseMatrix _global_M;
        SparseMatrix _global_G;
        std::vector<double> _global_b;

        std::vector<double> _solve;

        DomainFunctionsQuadsLinear _funcs;
        GridQuadLinear _grid;

        void generatePortrait();

        auto getLocalG(const MeshQuadLinear& mesh) -> LocalQuadsLinearMat;
        auto getLocalM(const MeshQuadLinear& mesh, bool addGamma = true) -> LocalQuadsLinearMat;
        auto getLocalB(const MeshQuadLinear& mesh) -> LocalQuadsLinearVec;
        
        void includeS2();
        void includeS1();

        auto findFinite(const Point& p) const -> size_t;
    };

}
