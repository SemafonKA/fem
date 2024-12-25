#include "three_dim_cuboid_linear_solver.h"

#include <format>

#include "../../logger.h"
#include "../../timer.h"
#include "../../gaussian_quadrature/gaussian_quadrature.h";
#include "../../three_steppers/Headers/IterSolvers.h";

using namespace std;

namespace fem::three_dim {

    template<size_t N>
    static inline std::array<double, N> operator* (double left, const std::array<double, N>& right) {
        std::array<double, N> res;
        for (size_t i = 0; i < N; i++) {
            res[i] = left * right[i];
        }
        return res;
    }

    static inline auto mu(size_t i) -> size_t {
        return i % 2;
    }

    static inline auto nu(size_t i) -> size_t {
        return (i / 2) % 2;
    }

    static inline auto v(size_t i) -> size_t {
        return i / 4;
    }

    static inline void addLocalMatrixToGlobal(const MeshCuboidLinear& rect, SparseMatrix& globalMat, const LocalCuboidLinearMat& localMat) {
        const auto& elems = rect.indOfPoints;

        for (int i = 0; i < elems.size(); i++) {
            // ��������� ��� ��������������� �������� �� ������ elems[i]
            for (int k = 0; k < i; k++) {
                // ���� ������� � ������� ��������������, �� �������
                if (elems[k] > elems[i]) {
                    continue;
                }

                auto id = globalMat.ig[elems[i]];
                while (id < globalMat.ig[elems[i] + 1] && globalMat.jg[id] != elems[k]) {
                    id++;
                }

                globalMat.ggl[id] += localMat[i][k];
                globalMat.ggu[id] += localMat[k][i];
            }
            // ��������� ������������ ��������
            globalMat.di[elems[i]] += localMat[i][i];
        }
    }

    static inline void addLocalbToGlobal(const MeshCuboidLinear& rect, std::vector<double>& globalVec, const LocalCuboidLinearVec& localVec) {
        const auto& elems = rect.indOfPoints;
        for (int i = 0; i < elems.size(); i++) {
            globalVec[elems[i]] += localVec[i];
        }
    }

    auto SolverCuboidLinear::value(Point p) -> double {
        auto mesh_ind = findFinite(p);
        if (mesh_ind == _grid.meshes.size()) {
            logger::error(format("Point ({}, {}, {}) is out of domain", p.x, p.y, p.z));
            return std::numeric_limits<double>::quiet_NaN();
        }
        const auto& mesh = _grid.meshes.at(mesh_ind);
        const array<Point, 2> pp{
            _grid.points[mesh.indOfPoints[0]],
            _grid.points[mesh.indOfPoints[7]],
        };

        auto basisLinear = [](size_t ind, double x, double x1, double x2) {
            if (ind == 0) {
                return (x2 - x) / (x2 - x1);
            }
            return (x - x1) / (x2 - x1);
            };
        auto psi = [&](size_t ind, double x, double y, double z) {
            auto xx = basisLinear(mu(ind), x, pp[0].x, pp[1].x);
            auto yy = basisLinear(nu(ind), y, pp[0].y, pp[1].y);
            auto zz = basisLinear(v(ind), z,  pp[0].z, pp[1].z);
            return xx * yy * zz;
            };

        double res = 0.0;
        for (size_t i = 0; i < 8; i++) {
            res += psi(i, p.x, p.y, p.z) * _solve.at(mesh.indOfPoints[i]);
        }
        return res;
    }

    [[nodiscard]]
    auto SolverCuboidLinear::solveStatic() -> std::vector<double> {
        Timer global_timer;
        Timer timer;
        global_timer.start();
        timer.start();
        generatePortrait();
        _global_G = SparseMatrix::copyShape(_global_mat);
        _global_M = SparseMatrix::copyShape(_global_mat);
        _global_b.resize(_global_mat.Size());
        timer.stop();
        logger::debug(format("SolverQuadsLinear::solveStatic - portrait generated by {} ms", timer.elapsedMilliseconds()));
        logger::debug(format("                                 size of portrait: {0} x {0}", _global_mat.Size()));

        timer.start();
        for (const auto& elem : _grid.meshes) {
            addLocalMatrixToGlobal(elem, _global_G, getLocalG(elem));
            addLocalMatrixToGlobal(elem, _global_M, getLocalM(elem));
            addLocalbToGlobal(elem, _global_b, getLocalB(elem));
        }
        timer.stop();
        logger::debug(format("SolverQuadsLinear::solveStatic - global matrices was computed by {} ms", timer.elapsedMilliseconds()));

        _global_mat = _global_M + _global_G;

        timer.start();
        includeS2();
        timer.stop();
        logger::debug(format("SolverQuadsLinear::solveStatic - s2 was computed by {} ms", timer.elapsedMilliseconds()));

        timer.start();
        includeS1();
        timer.stop();
        logger::debug(format("SolverQuadsLinear::solveStatic - s1 was computed by {} ms", timer.elapsedMilliseconds()));

        // TODO: ������ ��������� ���� ????

        timer.start();
        _solve.resize(_global_mat.Size());
        IterSolvers::LOS::Init_LuPrecond(_global_mat.Size(), _global_mat);
        IterSolvers::minEps = 1e-20;
        double eps;
        auto iter = IterSolvers::LOS::LuPrecond(_global_mat, _global_b, _solve, eps, false);
        IterSolvers::Destruct();
        timer.stop();
        logger::log(format("SolverQuadsLinear::solveStatic - SLAU solved by {} ms", timer.elapsedMilliseconds()));
        logger::log(format("                                 Iterations: {}, eps: {}", iter, eps));

        global_timer.stop();
        logger::log(format("SolverQuadsLinear::solveStatic - FEM solved by {} seconds", global_timer.elapsedSeconds()));
        return _solve;
    }

    void SolverCuboidLinear::generatePortrait() {
        const auto& nodes = _grid.points;

        _global_mat.di.resize(nodes.size());
        _global_mat.ig.resize(nodes.size() + 1);

        for (const auto& rect : _grid.meshes) {
            const auto& elems = rect.indOfPoints;
            for (int i = 0; i < elems.size(); i++) {
                for (int k = 0; k < i; k++) {
                    // ���� ������� � ������� ��������������, �� �������
                    if (elems[k] > elems[i])
                        continue;

                    bool isExist = false;
                    // ��������� �� ���� ������ ��� ��������, ���������� �� �����
                    // �������
                    for (auto it = _global_mat.ig[elems[i]]; it < _global_mat.ig[elems[i] + 1LL]; it++) {
                        if (_global_mat.jg[it] == elems[k]) {
                            isExist = true;
                            break;
                        }
                    }
                    if (!isExist) {
                        // ����, ���� �������� ������� ��������
                        auto it = _global_mat.ig[elems[i]];
                        while (it < _global_mat.ig[elems[i] + 1LL] && _global_mat.jg[it] < elems[k])
                            it++;

                        // ��� ������� ����� ����� �������� ������� �� ������, ���
                        // ���...
                        _global_mat.jg.insert(_global_mat.jg.begin() + it, elems[k]);

                        // ��������� ���� ��������� ig � ������� elems[i]+1 ����
                        // �������
                        for (auto j = elems[i] + 1; j < _global_mat.ig.size(); j++)
                            _global_mat.ig[j]++;
                    }
                }
            }
        }
        _global_mat.ggl.resize(_global_mat.jg.size());
        _global_mat.ggu.resize(_global_mat.jg.size());
    }

    auto SolverCuboidLinear::getLocalG(const MeshCuboidLinear& mesh) -> LocalCuboidLinearMat {
        LocalCuboidLinearMat g = {};
        constexpr std::array<double, 4> g1 = {
             1, -1,
            -1,  1
        };
        constexpr std::array<double, 4> m1 = {
            2.0 / 6.0, 1.0 / 6.0,
            1.0 / 6.0, 2.0 / 6.0,
        };
        const array<Point, 2> p{
            _grid.points[mesh.indOfPoints[0]],
            _grid.points[mesh.indOfPoints[7]],
        };
        auto lambda = _funcs.lambda(p[0].x, p[0].y, p[0].z, mesh.materialNum);

        double hx = p[1].x - p[0].x;
        double hy = p[1].y - p[0].y;
        double hz = p[1].z - p[0].z;

        auto gx = (1.0 / hx) * g1;
        auto gy = (1.0 / hy) * g1;
        auto gz = (1.0 / hz) * g1;

        auto mx = hx * m1;
        auto my = hy * m1;
        auto mz = hz * m1;

        for (size_t i = 0; i < 8; i++) {
            for (size_t j = 0; j < 8; j++) {
                double res = gx[mu(i) * 2 + mu(j)] * my[nu(i) * 2 + nu(j)] * mz[v(i) * 2 + v(j)];
                res += mx[mu(i) * 2 + mu(j)] * gy[nu(i) * 2 + nu(j)] * mz[v(i) * 2 + v(j)];
                res += mx[mu(i) * 2 + mu(j)] * my[nu(i) * 2 + nu(j)] * gz[v(i) * 2 + v(j)];
                res *= lambda;
                g[i][j] = res;
            }
        }

        return g;
    }

    auto SolverCuboidLinear::getLocalM(const MeshCuboidLinear& mesh, bool addGamma) -> LocalCuboidLinearMat {
        LocalCuboidLinearMat m = {};
        constexpr std::array<double, 4> m1 = {
            2.0 / 6.0, 1.0 / 6.0,
            1.0 / 6.0, 2.0 / 6.0,
        };
        const array<Point, 2> p{
            _grid.points[mesh.indOfPoints[0]],
            _grid.points[mesh.indOfPoints[7]],
        };
        auto gamma = addGamma ? _funcs.gamma(p[0].x, p[0].y, p[0].z, mesh.materialNum) : 1.0;

        double hx = p[1].x - p[0].x;
        double hy = p[1].y - p[0].y;
        double hz = p[1].z - p[0].z;

        auto mx = hx * m1;
        auto my = hy * m1;
        auto mz = hz * m1;

        for (size_t i = 0; i < 8; i++) {
            for (size_t j = 0; j < 8; j++) {
                double res = mx[mu(i) * 2 + mu(j)] * my[nu(i) * 2 + nu(j)] * mz[v(i) * 2 + v(j)];
                res *= gamma;
                m[i][j] = res;
            }
        }

        return m;
    }

    auto SolverCuboidLinear::getLocalB(const MeshCuboidLinear& mesh) -> LocalCuboidLinearVec {
        LocalCuboidLinearVec b = {};

        const array<Point, 8> p = {
            _grid.points[mesh.indOfPoints[0]],
            _grid.points[mesh.indOfPoints[1]],
            _grid.points[mesh.indOfPoints[2]],
            _grid.points[mesh.indOfPoints[3]],
            _grid.points[mesh.indOfPoints[4]],
            _grid.points[mesh.indOfPoints[5]],
            _grid.points[mesh.indOfPoints[6]],
            _grid.points[mesh.indOfPoints[7]],
        };

        auto m = getLocalM(mesh, false);
        for (size_t i = 0; i < 8; i++) {
            double sum = 0.0;
            for (size_t j = 0; j < 8; j++) {
                sum += m[i][j] * _funcs.func(p[j].x, p[j].y, p[j].z, mesh.materialNum);
            }
            b[i] = sum;
        }

        return b;
    }

    void SolverCuboidLinear::includeS2() {
        logger::warn("S2 boundaries are not implemented yet");
        // TODO: Implement
    }

    void SolverCuboidLinear::includeS1() {
        for (const auto& mesh : _grid.meshes) {
            const array<Point, 8> p = {
            _grid.points[mesh.indOfPoints[0]],
            _grid.points[mesh.indOfPoints[1]],
            _grid.points[mesh.indOfPoints[2]],
            _grid.points[mesh.indOfPoints[3]],
            _grid.points[mesh.indOfPoints[4]],
            _grid.points[mesh.indOfPoints[5]],
            _grid.points[mesh.indOfPoints[6]],
            _grid.points[mesh.indOfPoints[7]],
            };
            for (const auto& s1 : mesh.s1Conditions) {
                for (size_t i = 0; i < 4; i++) {
                    auto local_ind = s1.pointsIndices[i];
                    auto ind = mesh.indOfPoints[local_ind];
                    _global_mat.di.at(ind) = 1.0;
                    _global_b.at(ind) = _funcs.s1_func(p[local_ind].x, p[local_ind].y, p[local_ind].z, s1.functionNumber);

                    // Null string in the lower triangle part
                    for (auto j = _global_mat.ig[ind]; j < _global_mat.ig[ind + 1]; j++) {
                        _global_mat.ggl[j] = 0.0;
                    }

                    // Null string in the upper triangle part
                    for (size_t k = ind + 1; k < _global_mat.Size(); k++) {
                        for (size_t j = _global_mat.ig[k]; j < _global_mat.ig[k + 1]; j++) {
                            if (_global_mat.jg[j] == ind) {
                                _global_mat.ggu[j] = 0;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    auto SolverCuboidLinear::findFinite(const Point& p) const -> size_t {
        // Iterate through all meshes in the grid
        for (size_t i = 0; i < _grid.meshes.size(); ++i) {
            const auto& mesh = _grid.meshes[i];

            // Check if the point lies inside the mesh
            const array<Point, 2> pp{
                _grid.points[mesh.indOfPoints[0]],
                _grid.points[mesh.indOfPoints[7]],
            };

            bool isIn = true;
            isIn &= pp[0].x <= p.x;
            isIn &= pp[1].x >= p.x;
            isIn &= pp[0].y <= p.y;
            isIn &= pp[1].y >= p.y;
            isIn &= pp[0].z <= p.z;
            isIn &= pp[1].z >= p.z;
            if (isIn) {
                return i;
            }
        }

        return _grid.meshes.size(); // Return the size of the meshes vector if no mesh contains the point
    }
}