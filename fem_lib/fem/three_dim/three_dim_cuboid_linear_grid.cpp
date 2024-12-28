#include "three_dim_cuboid_linear_grid.h"

#include <vector>
#include <format>
#include <string>
#include <array>
#include <sstream>

#include "../../timer.h"
#include "../../logger.h"

using std::vector;
using std::array;
using std::format;
using std::string;
using fem::three_dim::Domain;
using std::ostringstream;


/**
 * @brief Returns Kx and Ky for grid (count of points by X and Y sides) regards to subdivisions and additional splits
 * @return [Kx, Ky]
 */
static auto getPointsCount(const Domain& domain) -> std::pair<size_t, size_t> {
    size_t Kx_ext = 1;
    size_t Ky_ext = 1;

    for (auto el : domain.nx) {
        Kx_ext += el * static_cast<size_t>(pow(2, domain.splitX));
    }
    for (auto el : domain.ny) {
        Ky_ext += el * static_cast<size_t>(pow(2, domain.splitY));
    }

    Kx_ext = Kx_ext == 1 ? 0 : Kx_ext;
    Ky_ext = Ky_ext == 1 ? 0 : Ky_ext;

    return { Kx_ext, Ky_ext };
}


/**
 * @brief Fill subdivide points for single axis and emplace it in grid by indices
 * @param grid - original grid needs to be changed (points must be preallocated)
 * @param beg - index of begin point of interval
 * @param end - index of end point of interval
 * @param step - step for index to get next point
 * @param cx - sparse coefficient
 * @param z - value of Z coordinate
 * @throw std::invalid_argument - in cases of: `beg == end`, `step == 0`, `cx <= 0`
 */
static void fillSubdivides(fem::three_dim::GridCuboidLinear& grid, size_t beg, size_t end, size_t step, double cx, double z) {
    // some debug checks
    if (beg == end) {
        auto error = format("Function `three_dim_quads_linear_grid.cpp::fillSubdivides`: parameters `beg` and `end` has equivalent values ({} and {})", beg, end);
        logger::debug(error);
        throw std::invalid_argument("Parameters `beg` and `end` must be different");
    }
    if (step == 0) {
        auto error = format("Function `three_dim_quads_linear_grid.cpp::fillSubdivides`: parameters `step` must be > 0, but it is {}.", step);
        logger::debug(error);
        throw std::invalid_argument("Parameter `step` must be > 0");
    }
    if (cx <= 0) {
        auto error = format("Function `three_dim_quads_linear_grid.cpp::fillSubdivides`: parameters `cx` must be > 0, but it is {}.", cx);
        logger::debug(error);
        throw std::invalid_argument("Parameter `cx` must be > 0");
    }

    auto xLenght = grid.points.at(end).x - grid.points.at(beg).x;
    auto yLenght = grid.points.at(end).y - grid.points.at(beg).y;
    size_t nx = (end - beg) / step;
    auto fraction = cx == 1.0 ? (1.0 / static_cast<double>(nx)) : ((cx - 1.0) / (std::pow(cx, nx) - 1));
    auto xh1 = xLenght * fraction;
    auto yh1 = yLenght * fraction;

    for (size_t i = 1; i < nx; i++) {
        const auto& p1 = grid.points.at(beg + (i - 1) * step);
        auto& p = grid.points.at(beg + i * step);

        p.x = p1.x + xh1 * std::pow(cx, i - 1);
        p.y = p1.y + yh1 * std::pow(cx, i - 1);
        p.z = z;
    }
}

/**
 * @brief Compute and fill point coordinates for chosen domain
 */
static void fillPoints(const Domain& domain, fem::three_dim::GridCuboidLinear& grid) {
    if (grid.Kx == 0 || grid.Ky == 0 || grid.Kz == 0) {
        auto warnMsg = format("Function `three_dim_quads_linear_grid.cpp::fillPoints`: some of sizes of grid (Kx, Ky, Kz) is equal to zero ({}, {}, {}). Maybe you need to set this values at first", grid.Kx, grid.Ky, grid.Kz);
        logger::debug(warnMsg, logger::Colors::warning);
        throw std::runtime_error("Bad grid size");
    }
    logger::log("XY-grid in domain will be rectangulated - only BOTTOM X values and only LEFT Y values will be used as point coordinates");

    auto z0 = domain.Z.at(0);
    size_t jj = 0; // Y index for grid
    for (size_t j = 0; j < domain.Ky - 1; j++) {
        size_t ii = 0; // X index for grid
        for (size_t i = 0; i < domain.Kx - 1; i++) {
            const array<size_t, 4> domainInd = {
                i + j * domain.Kx,           // Index of BL point
                i + 1 + j * domain.Kx,       // Index of BR point
                i + (j + 1) * domain.Kx,     // Index of TL point
                i + 1 + (j + 1) * domain.Kx, // Index of TR point
            };

            const array<size_t, 4> gridInd = {
                ii + jj * grid.Kx, // Index of BL point
                ii + (domain.nx[i] * static_cast<size_t>(pow(2, domain.splitX))) + jj * grid.Kx, // Index of BR point
                ii + (jj + domain.ny[j] * static_cast<size_t>(pow(2, domain.splitY))) * grid.Kx, // Index of TL point
                ii + (domain.nx[i] * static_cast<size_t>(pow(2, domain.splitX))) + (jj + domain.ny[j] * static_cast<size_t>(pow(2, domain.splitY))) * grid.Kx, // Index of TR point
            };

            // Set value for corner points
            for (size_t k = 0; k < gridInd.size(); k++) {
                auto& p = grid.points.at(gridInd[k]);
                p.x = domain.X.at(domainInd[k] % domain.Kx); // We got only BOTTOM X values
                p.y = domain.Y.at((domainInd[k] / domain.Kx) * domain.Kx); // and only LEFT Y values for rectangular grids
                p.z = z0;
            }

            // fill edges
            fillSubdivides(grid, gridInd[0], gridInd[1], 1, domain.cx[i], z0); // Bottom edge
            fillSubdivides(grid, gridInd[2], gridInd[3], 1, domain.cx[i], z0); // Top edge
            fillSubdivides(grid, gridInd[0], gridInd[2], grid.Kx, domain.cy[j], z0); // Left edge
            fillSubdivides(grid, gridInd[1], gridInd[3], grid.Kx, domain.cy[j], z0); // Right edge

            // fill inner points
            auto xStep = domain.nx[i] * static_cast<size_t>(pow(2, domain.splitX));
            for (size_t k = 1; k < xStep; k++) {
                fillSubdivides(grid, gridInd[0] + k, gridInd[2] + k, grid.Kx, domain.cy[j], z0);
            }

            // go to next ii
            ii += domain.nx[i] * static_cast<size_t>(pow(2, domain.splitX));
        }
        // go to next jj
        jj += domain.ny[j] * static_cast<size_t>(pow(2, domain.splitX));
    }

    // Extend nodes by Z axis
    const auto layerSize = grid.Kx * grid.Ky;
    for (size_t i = 1; i < grid.Kz; i++) {
        const auto layerBegin = i * layerSize; // Begin index of nodes by [i] Z axis
        for (size_t j = 0; j < layerSize; j++) {
            auto& gp = grid.points.at(j + layerBegin);
            gp = grid.points.at(j); // Copy nodes from [0] layer to [i] layer
            gp.z = domain.Z.at(i);
        }
    }
}

/**
 * @brief Generate and fill meshes for grid
 * @param grid - grid with pre-generated points
 */
static void fillMeshes(fem::three_dim::GridCuboidLinear& grid) {
    if (grid.Kx == 0 || grid.Ky == 0 || grid.Kz == 0) {
        auto warnMsg = format("Function `three_dim_quads_linear_grid.cpp::fillMeshes`: some of sizes of grid (Kx and Ky) is equal to zero ({} or {}). Maybe you need to use `three_dim_quads_linear_grid.cpp::fillPoints` at first", grid.Kx, grid.Ky);
        logger::debug(warnMsg, logger::Colors::warning);
    }

    const auto n = grid.Kx; // count of nodes by X
    const auto m = grid.Ky; // count of nodes by Y
    const auto k = grid.Kz; // count of nodes by Z 

    // Reserve memory for meshes in grid
    grid.meshes.reserve((n - 1) * (m - 1) * (k - 1));

    for (size_t z = 1; z < k; z++) {
        for (size_t y = 1; y < m; y++) {
            for (size_t x = 1; x < n; x++) {
                auto mesh = fem::three_dim::MeshCuboidLinear();
                mesh.indOfPoints = {
                    (x - 1) + (y - 1) * n + (z - 1) * n * m,
                    x + (y - 1) * n + (z - 1) * n * m,
                    (x - 1) + y * n + (z - 1) * n * m,
                    x + y * n + (z - 1) * n * m,
                    (x - 1) + (y - 1) * n + z * n * m,
                    x + (y - 1) * n + z * n * m,
                    (x - 1) + y * n + z * n * m,
                    x + y * n + z * n * m,
                };
                grid.meshes.push_back(std::move(mesh));
            }
        }
    }
}

/**
 * @brief Enumerate edges for every mesh element
 */
static void fillEdgeNumbers(fem::three_dim::GridCuboidLinear& grid) {
    throw std::runtime_error("Not implemented");
    // TODO: Implement it??

    auto ve = grid.Kx; // count of vertical edges
    auto he = ve - 1;  // count of horizontal edges
    auto elementsInRow = he;

    for (size_t i = 0; i < grid.meshes.size(); i++) {
        size_t row = i / elementsInRow; // current row of mesh elements
        auto& elem = grid.meshes.at(i);

        elem.indOfEdges[0] = elem.indOfPoints[0] + he * row;
        elem.indOfEdges[1] = elem.indOfEdges[0] + he;
        elem.indOfEdges[2] = elem.indOfEdges[1] + 1;
        elem.indOfEdges[3] = elem.indOfEdges[0] + he + ve;
    }
}

/**
 * @brief Fill materials of meshes by subdomains
 */
static void fillMaterials(const Domain& domain, fem::three_dim::GridCuboidLinear& grid) {
    if (grid.Kx == 0 || grid.Ky == 0 || grid.Kz == 0) {
        auto warnMsg = format("Function `three_dim_quads_linear_grid.cpp::fillMaterials`: some of sizes of grid (Kx and Ky) is equal to zero ({} or {}). Maybe you need to use `three_dim_quads_linear_grid.cpp::fillPoints` at first", grid.Kx, grid.Ky);
        logger::debug(warnMsg, logger::Colors::warning);
        return;
    }
    if (grid.meshes.empty()) {
        auto warnMsg = format("Function `three_dim_quads_linear_grid.cpp::fillMaterials`: there is no meshes in grid (size = {}). Maybe you need to use `three_dim_quads_linear_grid.cpp::fillMeshes` at first", grid.meshes.size());
        logger::debug(warnMsg, logger::Colors::warning);
        return;
    }

    // Domain to Grid points
    auto XlinesToPoints = std::vector<size_t>(domain.Kx);
    XlinesToPoints.at(0) = 0;
    for (size_t i = 0; i < domain.nx.size(); i++) {
        XlinesToPoints[i + 1] = XlinesToPoints[i] + domain.nx[i] * static_cast<size_t>(pow(2, domain.splitX));
    }
    auto YlinesToPoints = std::vector<size_t>(domain.Ky);
    YlinesToPoints.at(0) = 0;
    for (size_t i = 0; i < domain.ny.size(); i++) {
        YlinesToPoints[i + 1] = YlinesToPoints[i] + domain.ny[i] * static_cast<size_t>(pow(2, domain.splitY));
    }

    const auto kxky = grid.Kx * grid.Ky;

    // TODO: Refactor this to lower complexity from O(N2) to O(N logN)
    for (const auto& sub : domain.subdomains) {
        bool isSubInUse = false;
        // Convert subdomain coordinate lines to grid coordinate lines
        auto xBegLine = XlinesToPoints[sub.xBeginNum - 1];
        auto xEndLine = XlinesToPoints[sub.xEndNum - 1];
        auto yBegLine = YlinesToPoints[sub.yBeginNum - 1];
        auto yEndLine = YlinesToPoints[sub.yEndNum - 1];
        auto zBegLine = sub.zBeginNum - 1;
        auto zEndLine = sub.zEndNum - 1;

        for (auto& mesh : grid.meshes) {
            // Get coordinate lines for BL and TR nodes of mesh
            auto blLineX = mesh.indOfPoints[0] % kxky % grid.Kx;
            auto blLineY = mesh.indOfPoints[0] % kxky / grid.Kx;
            auto blLineZ = mesh.indOfPoints[0] / kxky;
            auto trLineX = mesh.indOfPoints[7] % kxky % grid.Kx;
            auto trLineY = mesh.indOfPoints[7] % kxky / grid.Kx;
            auto trLineZ = mesh.indOfPoints[7] / kxky;

            bool isMeshInSub = true;
            isMeshInSub &= blLineX >= xBegLine;
            isMeshInSub &= trLineX <= xEndLine;
            isMeshInSub &= blLineY >= yBegLine;
            isMeshInSub &= trLineY <= yEndLine;
            isMeshInSub &= blLineZ >= zBegLine;
            isMeshInSub &= trLineZ <= zEndLine;

            if (isMeshInSub) {
                mesh.materialNum = sub.materialNum;
                isSubInUse = true;
            }
        }
        if (isSubInUse) {
            grid.usedMaterials.push_back(sub.materialNum);
        }
        else {
            // TODO: Подумать, возможно ли такое вообще?
            // Область задаётся по координатным линиям. При этом координатные линии должны задаваться корректно. 
            // Единственное, что может пойти не так: некоторые области в итоге будут полностью перекрыты. 
            // Вот эти области в идеале надо исключать
        }
    }
}


/**
 * @brief Fill s1 boundary conditions for meshes
 */
static void fillS1(const Domain& domain, fem::three_dim::GridCuboidLinear& grid) {
    if (grid.Kx == 0 || grid.Ky == 0 || grid.Kz == 0) {
        auto warnMsg = format("Function `three_dim_quads_linear_grid.cpp::fillMaterials`: some of sizes of grid (Kx and Ky) is equal to zero ({} or {}). Maybe you need to use `three_dim_quads_linear_grid.cpp::fillPoints` at first", grid.Kx, grid.Ky);
        logger::debug(warnMsg, logger::Colors::warning);
        return;
    }
    if (grid.meshes.empty()) {
        auto warnMsg = format("Function `three_dim_quads_linear_grid.cpp::fillMaterials`: there is no meshes in grid (size = {}). Maybe you need to use `three_dim_quads_linear_grid.cpp::fillMeshes` at first", grid.meshes.size());
        logger::debug(warnMsg, logger::Colors::warning);
        return;
    }

    // Domain to Grid points
    auto XlinesToPoints = std::vector<size_t>(domain.Kx);
    XlinesToPoints.at(0) = 0;
    for (size_t i = 0; i < domain.nx.size(); i++) {
        XlinesToPoints[i + 1] = XlinesToPoints[i] + domain.nx[i] * static_cast<size_t>(pow(2, domain.splitX));
    }
    auto YlinesToPoints = std::vector<size_t>(domain.Ky);
    YlinesToPoints.at(0) = 0;
    for (size_t i = 0; i < domain.ny.size(); i++) {
        YlinesToPoints[i + 1] = YlinesToPoints[i] + domain.ny[i] * static_cast<size_t>(pow(2, domain.splitY));
    }

    auto conditions = std::vector<fem::three_dim::EdgeCondition>(6); // S1 edge conditions in order: Front, Back, Left, Right, Top, Bottom
    conditions[0].pointsIndices = std::array<size_t, 4> {0, 1, 4, 5};
    conditions[1].pointsIndices = std::array<size_t, 4> {2, 3, 6, 7};
    conditions[2].pointsIndices = std::array<size_t, 4> {0, 2, 4, 6};
    conditions[3].pointsIndices = std::array<size_t, 4> {1, 3, 5, 7};
    conditions[4].pointsIndices = std::array<size_t, 4> {4, 5, 6, 7};
    conditions[5].pointsIndices = std::array<size_t, 4> {0, 1, 2, 3};
    for (auto& cond : conditions) {
        cond.type = fem::three_dim::EdgeCondition::Type::First;
    }

    const auto kxky = grid.Kx * grid.Ky;
    // TODO: Refactor this to lower complexity from O(N2) to O(N logN)
    for (const auto& s1 : domain.s1_edges) {
        // Convert s1 coordinate lines to grid coordinate lines
        auto xBegLine = XlinesToPoints[s1.xBeginNum - 1];
        auto xEndLine = XlinesToPoints[s1.xEndNum - 1];
        auto yBegLine = YlinesToPoints[s1.yBeginNum - 1];
        auto yEndLine = YlinesToPoints[s1.yEndNum - 1];
        auto zBegLine = s1.zBeginNum - 1;
        auto zEndLine = s1.zEndNum - 1;

        // Fill function numbers
        for (auto& cond : conditions) {
            cond.functionNumber = s1.funcNum;
        }

        for (auto& mesh : grid.meshes) {
            // Get coordinate lines for BL and TR nodes of mesh
            auto blLineX = mesh.indOfPoints[0] % kxky % grid.Kx;
            auto blLineY = mesh.indOfPoints[0] % kxky / grid.Kx;
            auto blLineZ = mesh.indOfPoints[0] / kxky;
            auto trLineX = mesh.indOfPoints[7] % kxky % grid.Kx;
            auto trLineY = mesh.indOfPoints[7] % kxky / grid.Kx;
            auto trLineZ = mesh.indOfPoints[7] / kxky;

            // Front panel
            bool isInEdge = true;
            isInEdge &= blLineX >= xBegLine;
            isInEdge &= blLineX <= xEndLine;
            isInEdge &= blLineZ >= zBegLine;
            isInEdge &= blLineZ <= zEndLine;
            isInEdge &= blLineY == yBegLine && blLineY == yEndLine;
            if (isInEdge) {
                mesh.s1Conditions.push_back(conditions[0]);
            }

            // Back panel
            isInEdge = true;
            isInEdge &= blLineX >= xBegLine;
            isInEdge &= blLineX <= xEndLine;
            isInEdge &= blLineZ >= zBegLine;
            isInEdge &= blLineZ <= zEndLine;
            isInEdge &= trLineY == yBegLine && trLineY == yEndLine;
            if (isInEdge) {
                mesh.s1Conditions.push_back(conditions[1]);
            }

            // Left panel
            isInEdge = true;
            isInEdge &= blLineY >= yBegLine;
            isInEdge &= blLineY <= yEndLine;
            isInEdge &= blLineZ >= zBegLine;
            isInEdge &= blLineZ <= zEndLine;
            isInEdge &= blLineX == xBegLine && blLineX == xEndLine;
            if (isInEdge) {
                mesh.s1Conditions.push_back(conditions[2]);
            }

            // Right panel
            isInEdge = true;
            isInEdge &= blLineY >= yBegLine;
            isInEdge &= blLineY <= yEndLine;
            isInEdge &= blLineZ >= zBegLine;
            isInEdge &= blLineZ <= zEndLine;
            isInEdge &= trLineX == xBegLine && trLineX == xEndLine;
            if (isInEdge) {
                mesh.s1Conditions.push_back(conditions[3]);
            }

            // Top panel
            isInEdge = true;
            isInEdge &= blLineX >= xBegLine;
            isInEdge &= blLineX <= xEndLine;
            isInEdge &= blLineY >= yBegLine;
            isInEdge &= blLineY <= yEndLine;
            isInEdge &= trLineZ == zBegLine && trLineZ == zEndLine;
            if (isInEdge) {
                mesh.s1Conditions.push_back(conditions[4]);
            }

            // Bottom panel
            isInEdge = true;
            isInEdge &= blLineX >= xBegLine;
            isInEdge &= blLineX <= xEndLine;
            isInEdge &= blLineY >= yBegLine;
            isInEdge &= blLineY <= yEndLine;
            isInEdge &= blLineZ == zBegLine && blLineZ == zEndLine;
            if (isInEdge) {
                mesh.s1Conditions.push_back(conditions[5]);
            }
        }
    }
}


/**
 * @brief Fill s2 boundary conditions for meshes
 */
static void fillS2(const Domain& domain, fem::three_dim::GridCuboidLinear& grid) {
    if (grid.Kx == 0 || grid.Ky == 0) {
        auto warnMsg = format("Function `three_dim_quads_linear_grid.cpp::fillMaterials`: some of sizes of grid (Kx and Ky) is equal to zero ({} or {}). Maybe you need to use `three_dim_quads_linear_grid.cpp::fillPoints` at first", grid.Kx, grid.Ky);
        logger::debug(warnMsg, logger::Colors::warning);
    }
    if (grid.meshes.empty()) {
        auto warnMsg = format("Function `three_dim_quads_linear_grid.cpp::fillMaterials`: there is no meshes in grid (size = {}). Maybe you need to use `three_dim_quads_linear_grid.cpp::fillMeshes` at first", grid.meshes.size());
        logger::debug(warnMsg, logger::Colors::warning);
    }

    throw std::runtime_error("Not implemented");
    // TODO: Implement it ??

    // Domain to Grid points
    auto XlinesToPoints = std::vector<size_t>(domain.Kx);
    XlinesToPoints.at(0) = 0;
    for (size_t i = 0; i < domain.nx.size(); i++) {
        XlinesToPoints[i + 1] = XlinesToPoints[i] + domain.nx[i] * (domain.splitX + 1);
    }
    auto YlinesToPoints = std::vector<size_t>(domain.Ky);
    YlinesToPoints.at(0) = 0;
    for (size_t i = 0; i < domain.ny.size(); i++) {
        YlinesToPoints[i + 1] = YlinesToPoints[i] + domain.ny[i] * (domain.splitY + 1);
    }

    // TODO: Refactor this to lower complexity from O(N2) to O(N logN)
    for (const auto& s2 : domain.s2_edges) {
        // Convert s2 coordinate lines to grid coordinate lines
        auto xBegLine = XlinesToPoints[s2.xBeginNum - 1];
        auto xEndLine = XlinesToPoints[s2.xEndNum - 1];
        auto yBegLine = YlinesToPoints[s2.yBeginNum - 1];
        auto yEndLine = YlinesToPoints[s2.yEndNum - 1];

        // Determ if this edge vertical or horizontal
        auto isVert = xBegLine == xEndLine;
        auto s2_cond = fem::three_dim::EdgeCondition();
        s2_cond.type = fem::three_dim::EdgeCondition::Type::Second;
        s2_cond.functionNumber = s2.funcNum;

        for (auto& mesh : grid.meshes) {
            bool isHaveCondition = false;
            // Get coordinate lines for BL and TR nodes of mesh
            auto blLineX = mesh.indOfPoints[0] % grid.Kx;
            auto blLineY = mesh.indOfPoints[0] / grid.Kx;
            auto trLineX = mesh.indOfPoints[3] % grid.Kx;
            auto trLineY = mesh.indOfPoints[3] / grid.Kx;

            if (isVert) {
                auto isLeft = blLineX == xBegLine;
                auto isRight = trLineX == xEndLine;

                auto isInY = true;
                isInY &= blLineY >= yBegLine;
                isInY &= trLineY <= yEndLine;

                if (isInY && isLeft) {
                    s2_cond.pointsIndices[0] = 0;
                    s2_cond.pointsIndices[1] = 2;
                    isHaveCondition = true;
                }
                else if (isInY && isRight) {
                    s2_cond.pointsIndices[0] = 1;
                    s2_cond.pointsIndices[1] = 3;
                    isHaveCondition = true;
                }
            }
            else {
                auto isBottom = blLineY == yBegLine;
                auto isTop = trLineY == yEndLine;

                auto isInX = true;
                isInX &= blLineX >= xBegLine;
                isInX &= trLineX <= xEndLine;

                if (isInX && isBottom) {
                    s2_cond.pointsIndices[0] = 0;
                    s2_cond.pointsIndices[1] = 1;
                    isHaveCondition = true;
                }
                else if (isInX && isTop) {
                    s2_cond.pointsIndices[0] = 2;
                    s2_cond.pointsIndices[1] = 3;
                    isHaveCondition = true;
                }
            }

            if (isHaveCondition) {
                mesh.s2Conditions.push_back(s2_cond);
            }
        }
    }
}

/**
 * @brief Remove an unused meshes from grid (unused == doesn't contain any material)
 */
static void removeUnusedMeshes(fem::three_dim::GridCuboidLinear& grid) {
    auto& meshes = grid.meshes;

    auto filtered_meshes = std::vector<fem::three_dim::MeshCuboidLinear>();
    filtered_meshes.reserve(meshes.size());
    for (auto& mesh : meshes) {
        if (mesh.materialNum != 0) {
            filtered_meshes.push_back(std::move(mesh));
        }
    }

    grid.meshes = std::move(filtered_meshes);
}


namespace fem::three_dim {

    /**
     * @brief Build grid from provided domain
     * @param domain - the domain from which the grid will be built
     */
    void GridCuboidLinear::buildFrom(const Domain& domain) {
        auto timer = Timer(); // Debug timer

        // Get overall point count
        auto [countX, countY] = getPointsCount(domain);
        Kx = countX; Ky = countY; Kz = domain.Kz;
        points.resize(Kx * Ky * Kz);
        logger::debug(format("GridCuboidLinear::buildFrom: Grid will be {}x{}x{} size ({} overall)", Kx, Ky, Kz, points.size()));

        // Fill points
        timer.start();
        fillPoints(domain, *this);
        timer.stop();
        logger::debug(format("GridCuboidLinear::buildFrom: Points was computed in {} ms", timer.elapsedMilliseconds()));

        // Fill meshes
        timer.start();
        fillMeshes(*this);
        timer.stop();
        logger::debug(format("GridCuboidLinear::buildFrom: Meshes was computed in {} ms", timer.elapsedMilliseconds()));

        // Fill edges
        //timer.start();
        //fillEdgeNumbers(*this);
        //timer.stop();
        //logger::debug(format("GridCuboidLinear::buildFrom: Edges was enumerated in {} ms", timer.elapsedMilliseconds()));

        // Fill materials
        timer.start();
        fillMaterials(domain, *this);
        timer.stop();
        logger::debug(format("GridCuboidLinear::buildFrom: Materials was setted in {} ms", timer.elapsedMilliseconds()));

        // Fill s2
        //timer.start();
        //fillS2(domain, *this);
        //timer.stop();
        //logger::debug(format("GridCuboidLinear::buildFrom: S2 was setted in {} ms", timer.elapsedMilliseconds()));

        // Fill s1
        timer.start();
        fillS1(domain, *this);
        timer.stop();
        logger::debug(format("GridCuboidLinear::buildFrom: S1 was setted in {} ms", timer.elapsedMilliseconds()));

        // Remove unused meshes with empty materials
        timer.start();
        removeUnusedMeshes(*this);
        timer.stop();
        logger::debug(format("GridCuboidLinear::buildFrom: Unused meshes was removed in {} ms", timer.elapsedMilliseconds()));
    }

    /**
     * @brief Format grid points to string in format like `Kx Ky \\n x[0] y[0]   x[1] y[1]   ... \\n x[Kx] y[Kx]   ...`
     * @return formatted string
     */
    auto GridCuboidLinear::dumpPoints() const -> std::string {
        auto oss = ostringstream();

        oss << format("{} {}\n\n", Kx, Ky);

        for (size_t j = 0; j < Ky; j++) {
            for (size_t i = 0; i < Kx; i++) {
                const auto& pt = points[i + Kx * j];

                oss << format("{: 20.14e} {: 20.14e}   ", pt.x, pt.y);
            }
            oss << "\n";
        }

        return oss.str();
    }

    /**
     * @brief Format grid meshes to string in format like `Km \\n m[i]   p[0] p[1] p[2] p[3] \\n ...`
     * @return formatted string
     */
    auto GridCuboidLinear::dumpMeshes() const -> std::string {
        auto oss = ostringstream();

        oss << format("{}\n\n", meshes.size());

        for (const auto& mesh : meshes) {
            oss << format("{:3}   {:05} {:05} {:05} {:05}   {:05} {:05} {:05} {:05}\n",
                mesh.materialNum,
                mesh.indOfPoints[0], mesh.indOfPoints[1], mesh.indOfPoints[2], mesh.indOfPoints[3],
                mesh.indOfEdges[0], mesh.indOfEdges[1], mesh.indOfEdges[2], mesh.indOfEdges[3]);
        }

        return oss.str();
    }

    /**
     * @brief Dump grid to string if computer-readable format (see `dumpPoints` and `dumpMeshes` to know what format is)
     * @return formatted string
     */
    auto GridCuboidLinear::dump() const -> std::string {
        auto oss = ostringstream();
        oss << dumpPoints() << "\n";
        oss << dumpMeshes();
        return oss.str();
    }

}
