#include "two_dim_quads_linear_solver.h"

#include <format>

#include "../../logger.h"
#include "../../timer.h"
#include "../../gaussian_quadrature/gaussian_quadrature.h";
#include "../../three_steppers/Headers/IterSolvers.h";

using namespace std;


/**
 * @param xi - in range [0;1]
 * @param eta - in range [0;1]
 */
static inline auto x(double xi, double eta, const array<fem::two_dim::Point, 4>& points) -> double {
    double res = 0.0;
    res += (1 - xi) * (1 - eta) * points[0].x;
    res += xi * (1 - eta) * points[1].x;
    res += (1 - xi) * eta * points[2].x;
    res += xi * eta * points[3].x;
    return res;
}

/**
 * @param xi - in range [0;1]
 * @param eta - in range [0;1]
 */
static inline auto y(double xi, double eta, const array<fem::two_dim::Point, 4>& points) -> double {
    double res = 0.0;
    res += (1 - xi) * (1 - eta) * points[0].y;
    res += xi * (1 - eta) * points[1].y;
    res += (1 - xi) * eta * points[2].y;
    res += xi * eta * points[3].y;
    return res;
}


/**
 * @param val - xi or eta in range [0; 1]
 * @param ind - 0 for `1 - val` and 1 for `val`
 * @return
 */
static inline auto W(double val, size_t ind) -> double {
    if (ind == 0) {
        return 1 - val;
    }
    return val;
}

/**
 * @param val - xi or eta in range [0; 1]
 * @param ind - 0 for `-1` and 1 for `1`
 * @return
 */
static inline auto dW([[maybe_unused]] double val, size_t ind) -> double {
    if (ind == 0) {
        return -1;
    }
    return 1;
}

/**
 * @param xi - in range [0;1]
 * @param eta - in range [0;1]
 * @param ind - 0 for W1W1, 1 for W2W1, 2 for W1W2, 3 for W2W2
 */
static inline auto phi(double xi, double eta, size_t ind) -> double {
    auto xi_ind = ind % 2;
    auto eta_ind = ind / 2;

    return W(xi, xi_ind) * W(eta, eta_ind);
}

/**
 * @param xi - in range [0;1]
 * @param eta - in range [0;1]
 * @param ind - 0 for dW1W1, 1 for dW2W1, 2 for dW1W2, 3 for dW2W2
 */
static inline auto dPhi_dXi(double xi, double eta, size_t ind) -> double {
    auto xi_ind = ind % 2;
    auto eta_ind = ind / 2;

    return dW(xi, xi_ind) * W(eta, eta_ind);
}

/**
 * @param xi - in range [0;1]
 * @param eta - in range [0;1]
 * @param ind - 0 for W1dW1, 1 for W2dW1, 2 for W1dW2, 3 for W2dW2
 */
static inline auto dPhi_dEta(double xi, double eta, size_t ind) -> double {
    auto xi_ind = ind % 2;
    auto eta_ind = ind / 2;

    return W(xi, xi_ind) * dW(eta, eta_ind);
}

/**
 * @param xi - in range [0;1]
 * @param eta - in range [0;1]
 */
static inline auto J(double xi, double eta, const array<fem::two_dim::Point, 4>& points) -> double {
    std::array<double, 4> x = {
        points[0].x, points[1].x, points[2].x, points[3].x
    };
    std::array<double, 4> y = {
        points[0].y, points[1].y, points[2].y, points[3].y
    };

    double a0 = (x[1] - x[0]) * (y[2] - y[0]) - (y[1] - y[0]) * (x[2] - x[0]);
    double a1 = (x[1] - x[0]) * (y[3] - y[2]) - (y[1] - y[0]) * (x[3] - x[2]);
    double a2 = (y[2] - y[0]) * (x[3] - x[1]) - (x[2] - x[0]) * (y[3] - y[1]);

    return a0 + a1 * xi + a2 * eta;
}


namespace fem::two_dim {

    // Function to calculate the area of a triangle given its vertices
    static auto area(const Point& A, const Point& B, const Point& C) -> double {
        return 0.5 * abs((B.x - A.x) * (C.y - A.y) - (C.x - A.x) * (B.y - A.y));
    }

    // Function to check if a point is inside a triangle
    static auto isInTriangle(const Point& P, const Point& A, const Point& B, const Point& C) -> bool {
        double areaABC = area(A, B, C);
        double areaPAB = area(P, A, B);
        double areaPBC = area(P, B, C);
        double areaPCA = area(P, C, A);

        // If the sum of the areas of sub-triangles equals the area of the main triangle,
        // then the point is inside the triangle.
        return (abs(areaABC - (areaPAB + areaPBC + areaPCA)) < 1e-9);
    }

    static inline void addLocalMatrixToGlobal(const MeshQuadLinear& rect, SparseMatrix& globalMat, const LocalQuadsLinearMat& localMat) {
        const auto& elems = rect.indOfPoints;

        for (int i = 0; i < 4; i++) {
            // добавляем все внедиагональные элементы на строке elems[i]
            for (int k = 0; k < i; k++) {
                // Если элемент в верхнем прямоугольнике, то скипаем
                if (elems[k] > elems[i]) {
                    continue;
                }

                auto id = globalMat.ig[elems[i]];
                while (id < globalMat.ig[elems[i] + 1ll] && globalMat.jg[id] != elems[k]) {
                    id++;
                }

                globalMat.ggl[id] += localMat[i][k];
                globalMat.ggu[id] += localMat[k][i];
            }
            // добавляем диагональные элементы
            globalMat.di[elems[i]] += localMat[i][i];
        }
    }

    static inline void addLocalbToGlobal(const MeshQuadLinear& rect, std::vector<double>& globalVec, const LocalQuadsLinearVec& localVec) {
        const auto& elems = rect.indOfPoints;
        for (int i = 0; i < 4; i++) {
            globalVec[elems[i]] += localVec[i];
        }
    }

    auto SolverQuadsLinear::value(Point p) -> double {
        auto mesh_ind = findFinite(p);
        if (mesh_ind == _grid.meshes.size()) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        const auto& mesh = _grid.meshes.at(mesh_ind);

        const array<Point, 4> mp = {
            _grid.points[mesh.indOfPoints[0]],
            _grid.points[mesh.indOfPoints[1]],
            _grid.points[mesh.indOfPoints[2]],
            _grid.points[mesh.indOfPoints[3]],
        };
        std::array<double, 4> x = {
            mp[0].x, mp[1].x, mp[2].x, mp[3].x
        };
        std::array<double, 4> y = {
            mp[0].y, mp[1].y, mp[2].y, mp[3].y
        };

        const array<double, 3> alpha = {
            (x[1] - x[0])* (y[2] - y[0]) - (y[1] - y[0]) * (x[2] - x[0]),
            (x[1] - x[0])* (y[3] - y[2]) - (y[1] - y[0]) * (x[3] - x[2]),
            (y[2] - y[0])* (x[3] - x[1]) - (x[2] - x[0]) * (y[3] - y[1])
        };
        const array<double, 7> beta = {
            0,
            x[2] - x[0],
            x[1] - x[0],
            y[2] - y[0],
            y[1] - y[0],
            x[0] - x[1] - x[2] + x[3],
            y[0] - y[1] - y[2] + y[3],
        };
        
        double w = beta[6] * (p.x - x[0]) - beta[5] * (p.y - y[0]);

        double xi = 0.0;
        double eta = 0.0;

        if (abs(alpha[1]) < 1e-10 && abs(alpha[2]) < 1e-10) {
            xi = beta[3] * (p.x - x[0]) - beta[1] * (p.y - y[0]);
            xi /= beta[2] * beta[3] - beta[1] * beta[4];

            eta = beta[2] * (p.y - y[0]) - beta[4] * (p.x - x[0]);
            eta /= beta[2] * beta[3] - beta[1] * beta[4];
        }
        else if (abs(alpha[1]) < 1e-10) {
            xi = alpha[2] * (p.x - x[0]) + beta[1] * w;
            xi /= alpha[2] * beta[2] - beta[5] * w;

            eta = -w / alpha[2];
        }
        else if (abs(alpha[2]) < 1e-10) {
            xi = w / alpha[1];

            eta = alpha[1] * (p.y - y[0]) - beta[4] * w;
            eta /= alpha[1] * beta[3] + beta[6] * w;
        }
        else {
            double a = beta[5] * alpha[2];
            double b = alpha[2] * beta[2] + alpha[1] * beta[1] + beta[5] * w;
            double c = alpha[1] * (x[0] - p.x) + beta[2] * w;

            double discriminant = b * b - 4 * a * c;
            if (discriminant < 0) {
                return std::numeric_limits<double>::quiet_NaN(); // No real roots
            }

            // Calculate the two possible solutions
            double sqrtDiscriminant = std::sqrt(discriminant);
            double eta1 = (-b + sqrtDiscriminant) / (2 * a);
            double eta2 = (-b - sqrtDiscriminant) / (2 * a);
            double xi1 = alpha[2] / alpha[1] * eta1 + w / alpha[1];
            double xi2 = alpha[2] / alpha[1] * eta2 + w / alpha[1];

            if (0 <= eta1 && eta1 <= 1 && 0 <= xi1 && xi1 <= 1) {
                xi = xi1;
                eta = eta1;
            }
            else {
                xi = xi2;
                eta = eta2;
            }
        }

        double res = 0.0;
        for (size_t i = 0; i < 4; i++) {
            res += phi(xi, eta, i) * _solve.at(mesh.indOfPoints[i]);
        }
        return res;
    }

    [[nodiscard]]
    auto SolverQuadsLinear::solveStatic() -> std::vector<double> {
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

        // TODO: Учесть фиктивные узлы ????

        timer.start();
        _solve.resize(_global_mat.Size());
        IterSolvers::LOS::Init_LuPrecond(_global_mat.Size(), _global_mat);
        IterSolvers::minEps = 1e-40;
        double eps;
        auto iter = IterSolvers::LOS::LuPrecond(_global_mat, _global_b, _solve, eps, false);
        IterSolvers::Destruct();
        timer.stop();
        logger::debug(format("SolverQuadsLinear::solveStatic - SLAU solved by {} ms", timer.elapsedMilliseconds()));
        logger::debug(format("                                 Iterations: {}, eps: {}", iter, eps));

        global_timer.stop();
        logger::log(format("SolverQuadsLinear::solveStatic - FEM solved by {} second", global_timer.elapsedSeconds()));
        return _solve;
    }

    void SolverQuadsLinear::generatePortrait() {
        const auto& nodes = _grid.points;

        _global_mat.di.resize(nodes.size());
        _global_mat.ig.resize(nodes.size() + 1);

        for (const auto& rect : _grid.meshes) {
            const auto& elems = rect.indOfPoints;
            for (int i = 0; i < elems.size(); i++) {
                for (int k = 0; k < i; k++) {
                    // Если элемент в верхнем прямоугольнике, то скипаем
                    if (elems[k] > elems[i])
                        continue;

                    bool isExist = false;
                    // Пробегаем по всей строке для проверки, существует ли такой
                    // элемент
                    for (auto it = _global_mat.ig[elems[i]]; it < _global_mat.ig[elems[i] + 1LL]; it++) {
                        if (_global_mat.jg[it] == elems[k]) {
                            isExist = true;
                            break;
                        }
                    }
                    if (!isExist) {
                        // Ищем, куда вставить элемент портрета
                        auto it = _global_mat.ig[elems[i]];
                        while (it < _global_mat.ig[elems[i] + 1LL] && _global_mat.jg[it] < elems[k])
                            it++;

                        // Для вставки нужно взять итератор массива от начала, так
                        // что...
                        _global_mat.jg.insert(_global_mat.jg.begin() + it, elems[k]);

                        // Добавляем всем элементам ig с позиции elems[i]+1 один
                        // элемент
                        for (auto j = elems[i] + 1; j < _global_mat.ig.size(); j++)
                            _global_mat.ig[j]++;
                    }
                }
            }
        }
        _global_mat.ggl.resize(_global_mat.jg.size());
        _global_mat.ggu.resize(_global_mat.jg.size());
    }

    auto SolverQuadsLinear::getLocalG(const MeshQuadLinear& mesh) -> LocalQuadsLinearMat {
        LocalQuadsLinearMat g = {};

        const array<Point, 4> p = {
            _grid.points[mesh.indOfPoints[0]],
            _grid.points[mesh.indOfPoints[1]],
            _grid.points[mesh.indOfPoints[2]],
            _grid.points[mesh.indOfPoints[3]],
        };

        const array<double, 6> beta = {
            p[2].x - p[0].x,
            p[1].x - p[0].x,
            p[2].y - p[0].y,
            p[1].y - p[0].y,
            p[0].x - p[1].x - p[2].x + p[3].x,
            p[0].y - p[1].y - p[2].y + p[3].y,
        };

        size_t i = 0;
        size_t j = 0;

        auto int1_func = [&](double xi, double eta) {
            auto _x = x(xi, eta, p);
            auto _y = y(xi, eta, p);
            double result = _funcs.lambda(_x, _y, mesh.materialNum);
            result *= dPhi_dXi(xi, eta, i) * (beta[5] * xi + beta[2]) - dPhi_dEta(xi, eta, i) * (beta[5] * eta + beta[3]);
            result *= dPhi_dXi(xi, eta, j) * (beta[5] * xi + beta[2]) - dPhi_dEta(xi, eta, j) * (beta[5] * eta + beta[3]);
            result /= std::abs(J(xi, eta, p));
            return result;
            };

        auto int2_func = [&](double xi, double eta) {
            auto _x = x(xi, eta, p);
            auto _y = y(xi, eta, p);
            double result = _funcs.lambda(_x, _y, mesh.materialNum);
            result *= dPhi_dEta(xi, eta, i) * (beta[4] * eta + beta[1]) - dPhi_dXi(xi, eta, i) * (beta[4] * xi + beta[0]);
            result *= dPhi_dEta(xi, eta, j) * (beta[4] * eta + beta[1]) - dPhi_dXi(xi, eta, j) * (beta[4] * xi + beta[0]);
            result /= std::abs(J(xi, eta, p));
            return result;
            };

        auto int1_solver = Gaussian_4p::TwoDimentionalSolver(0, 1, 0, 1, int1_func);
        auto int2_solver = Gaussian_4p::TwoDimentionalSolver(0, 1, 0, 1, int2_func);

        for (i = 0; i < 4; i++) {
            for (j = 0; j < 4; j++) {
                g[i][j] = int1_solver.compute() + int2_solver.compute();   // [i] & [j] variables are linked to functions
            }
        }

        return g;
    }

    auto SolverQuadsLinear::getLocalM(const MeshQuadLinear& mesh, bool addGamma) -> LocalQuadsLinearMat {
        LocalQuadsLinearMat m = {};

        const array<Point, 4> p = {
            _grid.points[mesh.indOfPoints[0]],
            _grid.points[mesh.indOfPoints[1]],
            _grid.points[mesh.indOfPoints[2]],
            _grid.points[mesh.indOfPoints[3]],
        };

        size_t i = 0;
        size_t j = 0;

        auto int_func = [&](double xi, double eta) {
            auto _x = x(xi, eta, p);
            auto _y = y(xi, eta, p);
            double result = addGamma ? _funcs.gamma(_x, _y, mesh.materialNum) : 1;
            result *= phi(xi, eta, i) * phi(xi, eta, j);
            result *= std::abs(J(xi, eta, p));
            return result;
            };

        auto int_solver = Gaussian_4p::TwoDimentionalSolver(0, 1, 0, 1, int_func);

        for (i = 0; i < 4; i++) {
            for (j = 0; j < 4; j++) {
                m[i][j] = int_solver.compute();   // [i] & [j] variables are linked to [int_func] function
            }
        }

        return m;
    }

    auto SolverQuadsLinear::getLocalB(const MeshQuadLinear& mesh) -> LocalQuadsLinearVec {
        LocalQuadsLinearVec b = {};

        const array<Point, 4> p = {
            _grid.points[mesh.indOfPoints[0]],
            _grid.points[mesh.indOfPoints[1]],
            _grid.points[mesh.indOfPoints[2]],
            _grid.points[mesh.indOfPoints[3]],
        };

        size_t i = 0;

        auto int_func = [&](double xi, double eta) {
            double result = 1;
            auto _x = x(xi, eta, p);
            auto _y = y(xi, eta, p);
            result *= phi(xi, eta, i);
            result *= std::abs(J(xi, eta, p)) * _funcs.func(_x, _y, mesh.materialNum);
            return result;
            };

        auto int_solver = Gaussian_4p::TwoDimentionalSolver(0, 1, 0, 1, int_func);

        for (i = 0; i < 4; i++) {
            b[i] = int_solver.compute();   // [i] variable are linked to [int_solver] function
        }

        // if ^ not work:
        //auto m = getLocalM(mesh, false);
        //for (size_t i = 0; i < 4; i++) {
        //    double sum = 0.0;
        //    for (size_t j = 0; j < 4; j++) {
        //        sum += m[i][j] * _funcs.func(p[j].x, p[j].y, mesh.materialNum);
        //    }
        //    b[i] = sum;
        //}

        return b;
    }

    void SolverQuadsLinear::includeS2() {
        for (const auto& mesh : _grid.meshes) {
            const array<Point, 4> p = {
                _grid.points[mesh.indOfPoints[0]],
                _grid.points[mesh.indOfPoints[1]],
                _grid.points[mesh.indOfPoints[2]],
                _grid.points[mesh.indOfPoints[3]],
            };

            for (const auto& s2 : mesh.s2Conditions) {
                // if edge is vertical
                if (s2.pointsIndices[0] == 0 && s2.pointsIndices[1] == 2 ||
                    s2.pointsIndices[0] == 1 && s2.pointsIndices[1] == 3) {
                    auto xi = s2.pointsIndices[0] == 0 ? 0.0 : 1.0;

                    size_t i = 0;

                    auto int_func = [&](double eta) {
                        double result = 1;
                        auto _x = x(xi, eta, p);
                        auto _y = y(xi, eta, p);
                        result *= phi(xi, eta, i);
                        result *= std::abs(J(xi, eta, p)) * _funcs.s2_func(_x, _y, s2.functionNumber);
                        return result;
                        };
                    auto int_solver = Gaussian_4p::OneDimentionSolver(0, 1, int_func);

                    for (size_t j = 0; j < 2; j++) {
                        i = s2.pointsIndices[j];
                        _global_b.at(mesh.indOfPoints[i]) += int_solver.compute();
                    }
                }
                // if edge is horizontal
                else {
                    auto eta = s2.pointsIndices[0] == 0 ? 0.0 : 1.0;

                    size_t i = 0;

                    auto int_func = [&](double xi) {
                        double result = 1;
                        auto _x = x(xi, eta, p);
                        auto _y = y(xi, eta, p);
                        result *= phi(xi, eta, i);
                        result *= std::abs(J(xi, eta, p)) * _funcs.s2_func(_x, _y, s2.functionNumber);
                        return result;
                        };
                    auto int_solver = Gaussian_4p::OneDimentionSolver(0, 1, int_func);

                    for (size_t j = 0; j < 2; j++) {
                        i = s2.pointsIndices[j];
                        _global_b.at(mesh.indOfPoints[i]) += int_solver.compute();
                    }
                }
            }
        }
    }

    void SolverQuadsLinear::includeS1() {
        for (const auto& mesh : _grid.meshes) {
            const array<Point, 4> p = {
                _grid.points[mesh.indOfPoints[0]],
                _grid.points[mesh.indOfPoints[1]],
                _grid.points[mesh.indOfPoints[2]],
                _grid.points[mesh.indOfPoints[3]],
            };

            for (const auto& s1 : mesh.s1Conditions) {
                for (size_t i = 0; i < 2; i++) {
                    auto local_ind = s1.pointsIndices[i];
                    auto ind = mesh.indOfPoints[local_ind];
                    _global_mat.di.at(ind) = 1.0;
                    _global_b.at(ind) = _funcs.s1_func(p[local_ind].x, p[local_ind].y, s1.functionNumber);

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

    auto SolverQuadsLinear::findFinite(const Point& p) const -> size_t {
        // Iterate through all meshes in the grid
        for (size_t i = 0; i < _grid.meshes.size(); ++i) {
            const auto& mesh = _grid.meshes[i];

            // Check if the point lies inside the mesh
            Point p1 = _grid.points[mesh.indOfPoints[0]];
            Point p2 = _grid.points[mesh.indOfPoints[1]];
            Point p3 = _grid.points[mesh.indOfPoints[2]];
            Point p4 = _grid.points[mesh.indOfPoints[3]];

            if (isInTriangle(p, p1, p2, p3) || isInTriangle(p, p1, p3, p4)) {
                return i; // Return the index of the mesh
            }
        }

        return _grid.meshes.size(); // Return the size of the meshes vector if no mesh contains the point
    }
}
