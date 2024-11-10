#include "two_dim_quads_linear_grid.h"

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
using fem::two_dim::Domain;
using std::ostringstream;


/**
 * @brief Returns Kx and Ky for grid (count of points by X and Y sides) regards to subdivisions and additional splits
 * @return [Kx, Ky]
 */
static auto getPointsCount(const Domain& domain) -> std::pair<size_t, size_t> {
    size_t Kx_ext = domain.Kx;
    size_t Ky_ext = domain.Ky;

    for (auto el : domain.nx) {
        Kx_ext += el - 1;
    }
    for (auto el : domain.ny) {
        Ky_ext += el - 1;
    }

    Kx_ext = Kx_ext * (domain.splitX + 1);
    Ky_ext = Ky_ext * (domain.splitX + 1);

    return { Kx_ext, Ky_ext };
}


/**
 * @brief Fill subdivide points for single axis and emplace it in grid by indices
 * @param grid - original grid needs to be changed (points must be preallocated)
 * @param beg - index of begin point of interval
 * @param end - index of end point of interval
 * @param step - step for index to get next point
 * @param cx - sparse coefficient
 * @throw std::invalid_argument - in cases of: `beg == end`, `step == 0`, `cx <= 0`
 */
static void fillSubdivides(fem::two_dim::GridQuadLinear& grid, size_t beg, size_t end, size_t step, double cx) {
    // some debug checks
    if (beg == end) {
        auto error = format("Function `two_dim_quads_linear_grid.cpp::fillSubdivides`: parameters `beg` and `end` has equivalent values ({} and {})", beg, end);
        logger::debug(error);
        throw std::invalid_argument("Parameters `beg` and `end` must be different");
    }
    if (step == 0) {
        auto error = format("Function `two_dim_quads_linear_grid.cpp::fillSubdivides`: parameters `step` must be > 0, but it is {}.", step);
        logger::debug(error);
        throw std::invalid_argument("Parameter `step` must be > 0");
    }
    if (cx <= 0) {
        auto error = format("Function `two_dim_quads_linear_grid.cpp::fillSubdivides`: parameters `cx` must be > 0, but it is {}.", cx);
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
        grid.points.at(beg + i * step).x = grid.points.at(beg + (i - 1) * step).x + xh1 * std::pow(cx, i - 1);
        grid.points.at(beg + i * step).y = grid.points.at(beg + (i - 1) * step).y + yh1 * std::pow(cx, i - 1);
    }
}

/**
 * @brief Compute and fill point coordinates for chosen domain
 */
static void fillPoints(const Domain& domain, fem::two_dim::GridQuadLinear& grid) {
    if (grid.Kx == 0 || grid.Ky == 0) {
        auto warnMsg = format("Function `two_dim_quads_linear_grid.cpp::fillPoints`: some of sizes of grid (Kx and Ky) is equal to zero ({} or {}). Maybe you need to use `two_dim_quads_linear_grid.cpp::fillPoints` at first", grid.Kx, grid.Ky);
        logger::debug(warnMsg, logger::Colors::warning);
        throw std::runtime_error("Bad grid size");
    }

    size_t jj = 0; // Y index for grid
    for (size_t j = 0; j < domain.Ky - 1; j++) {
        size_t ii = 0; // X index for grid
        for (size_t i = 0; i < domain.Kx - 1; i++) {
            array<size_t, 4> domainInd = {
                i + j * domain.Kx,           // Index of BL point
                i + 1 + j * domain.Kx,       // Index of BR point
                i + (j + 1) * domain.Kx,     // Index of TL point
                i + 1 + (j + 1) * domain.Kx, // Index of TR point
            };

            array<size_t, 4> gridInd = {
                ii + jj * grid.Kx, // Index of BL point
                ii + (domain.nx[i] * (domain.splitX + 1)) + jj * grid.Kx, // Index of BR point
                ii + (jj + domain.ny[j] * (domain.splitY + 1)) * grid.Kx, // Index of TL point
                ii + (domain.nx[i] * (domain.splitX + 1)) + (jj + domain.ny[j] * (domain.splitY + 1)) * grid.Kx, // Index of TR point
            };

            // Set value for corner points
            for (size_t k = 0; k < gridInd.size(); k++) {
                grid.points.at(gridInd[k]).x = domain.X.at(domainInd[k]);
                grid.points.at(gridInd[k]).y = domain.Y.at(domainInd[k]);
            }

            // fill edges
            fillSubdivides(grid, gridInd[0], gridInd[1], 1, domain.cx[i]); // Bottom edge
            fillSubdivides(grid, gridInd[2], gridInd[3], 1, domain.cx[i]); // Top edge
            fillSubdivides(grid, gridInd[0], gridInd[2], grid.Kx, domain.cy[j]); // Left edge
            fillSubdivides(grid, gridInd[1], gridInd[3], grid.Kx, domain.cy[j]); // Right edge

            // fill inner points
            auto xStep = domain.nx[i] * (domain.splitX + 1);
            for (size_t k = 1; k < xStep; k++) {
                fillSubdivides(grid, gridInd[0] + k, gridInd[2] + k, grid.Kx, domain.cy[j]);
            }

            // go to next ii
            ii += domain.nx[i] * (domain.splitX + 1);
        }
        // go to next jj
        jj += domain.ny[j] * (domain.splitY + 1);
    }
}

/**
 * @brief Generate and fill meshes for grid
 * @param grid - grid with pre-generated points
 */
static void fillMeshes(fem::two_dim::GridQuadLinear& grid) {
    if (grid.Kx == 0 || grid.Ky == 0) {
        auto warnMsg = format("Function `two_dim_quads_linear_grid.cpp::fillMeshes`: some of sizes of grid (Kx and Ky) is equal to zero ({} or {}). Maybe you need to use `two_dim_quads_linear_grid.cpp::fillPoints` at first", grid.Kx, grid.Ky);
        logger::debug(warnMsg, logger::Colors::warning);
    }

    auto n = grid.Kx; // count of nodes by X
    auto m = grid.Ky; // count of nodes by Y

    // Reserve memory for meshes in grid
    grid.meshes.reserve((n - 1) * (m - 1));

    for (size_t j = 0; j < m - 1; j++) {
        for (size_t i = 0; i < n - 1; i++) {
            auto mesh = fem::two_dim::MeshQuadLinear();
            mesh.indOfPoints = {
                i + j * n,
                (i + 1) + j * n,
                i + (j + 1) * n,
                (i + 1) + (j + 1) * n,
            };
            grid.meshes.push_back(std::move(mesh));
        }
    }
}

/**
 * @brief Fill materials of meshes by subdomains
 */
static void fillMaterials(const Domain& domain, fem::two_dim::GridQuadLinear& grid) {
    if (grid.Kx == 0 || grid.Ky == 0) {
        auto warnMsg = format("Function `two_dim_quads_linear_grid.cpp::fillMaterials`: some of sizes of grid (Kx and Ky) is equal to zero ({} or {}). Maybe you need to use `two_dim_quads_linear_grid.cpp::fillPoints` at first", grid.Kx, grid.Ky);
        logger::debug(warnMsg, logger::Colors::warning);
    }
    if (grid.meshes.empty()) {
        auto warnMsg = format("Function `two_dim_quads_linear_grid.cpp::fillMaterials`: there is no meshes in grid (size = {}). Maybe you need to use `two_dim_quads_linear_grid.cpp::fillMeshes` at first", grid.meshes.size());
        logger::debug(warnMsg, logger::Colors::warning);
    }

    // Domain to Grid points
    auto XlinesToPoints = std::vector<size_t>(domain.Kx);
    XlinesToPoints.at(0) = 0;
    for (size_t i = 0; i < domain.nx.size(); i++) {
        XlinesToPoints[i + 1] = XlinesToPoints[i] + domain.nx[i];
    }
    auto YlinesToPoints = std::vector<size_t>(domain.Ky);
    YlinesToPoints.at(0) = 0;
    for (size_t i = 0; i < domain.ny.size(); i++) {
        YlinesToPoints[i + 1] = YlinesToPoints[i] + domain.ny[i];
    }

    for (const auto& sub : domain.subdomains) {
        bool isSubInUse = false;
        // Convert subdomain coordinate lines to grid coordinate lines
        auto xBegLine = XlinesToPoints[sub.xBeginNum - 1];
        auto xEndLine = XlinesToPoints[sub.xEndNum - 1];
        auto yBegLine = YlinesToPoints[sub.yBeginNum - 1];
        auto yEndLine = YlinesToPoints[sub.yEndNum - 1];

        for (auto& mesh : grid.meshes) {
            // Get coordinate lines for BL and TR nodes of mesh
            auto blLineX = mesh.indOfPoints[0] % grid.Kx;
            auto blLineY = mesh.indOfPoints[0] / grid.Kx;
            auto trLineX = mesh.indOfPoints[3] % grid.Kx;
            auto trLineY = mesh.indOfPoints[3] / grid.Kx;

            bool isMeshInSub = true;
            isMeshInSub &= blLineX >= xBegLine;
            isMeshInSub &= trLineX <= xEndLine;
            isMeshInSub &= blLineY >= yBegLine;
            isMeshInSub &= trLineY <= yEndLine;

            if (isMeshInSub) {
                mesh.materialNum = sub.materialNum;
                isSubInUse = true;
            }
        }
        if (isSubInUse) {
            grid.usedMaterials.push_back(sub.materialNum);
        }
    }

    // Additional checks
    bool allHasMaterial = true;
    for (const auto& mesh : grid.meshes) {
        allHasMaterial &= mesh.materialNum != 0;
    }
    if (!allHasMaterial) {
        logger::log("Some meshes has not material. Maybe you make some mistake with subdomains?", logger::Colors::warning);
    }
}


namespace fem::two_dim {

    /**
     * @brief Build grid from provided domain
     * @param domain - the domain from which the grid will be built
     */
    void GridQuadLinear::buildFrom(const Domain& domain) {
        auto timer = Timer(); // Debug timer

        // TODO: Throw bad_domain exception if domain is bad

        // Get overall point count
        auto [countX, countY] = getPointsCount(domain);
        Kx = countX; Ky = countY;
        points.resize(Kx * Ky);
        logger::debug(format("GridQuadLinear::buildFrom: Grid will be {}x{} size ({} overall)", Kx, Ky, points.size()));

        // Fill points
        timer.start();
        fillPoints(domain, *this);
        timer.stop();
        logger::debug(format("GridQuadLinear::buildFrom: Points was computed in {} ms", timer.elapsedMilliseconds()));

        // Fill meshes
        timer.start();
        fillMeshes(*this);
        timer.stop();
        logger::debug(format("GridQuadLinear::buildFrom: Meshes was computed in {} ms", timer.elapsedMilliseconds()));

        // Fill materials
        timer.start();
        fillMaterials(domain, *this);
        timer.stop();
        logger::debug(format("GridQuadLinear::buildFrom: Materials was setted in {} ms", timer.elapsedMilliseconds()));
    }

    /**
     * @brief Format grid points to string in format like `Kx Ky \\n x[0] y[0]   x[1] y[1]   ... \\n x[Kx] y[Kx]   ...`
     * @return formatted string
     */
    auto GridQuadLinear::dumpPoints() const -> std::string {
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
    auto GridQuadLinear::dumpMeshes() const -> std::string {
        auto oss = ostringstream();

        oss << format("{}\n\n", meshes.size());

        for (const auto& mesh : meshes) {
            oss << format("{:3}   {:05} {:05} {:05} {:05}\n", mesh.materialNum, mesh.indOfPoints[0], mesh.indOfPoints[1], mesh.indOfPoints[2], mesh.indOfPoints[3]);
        }

        return oss.str();
    }

    /**
     * @brief Dump grid to string if computer-readable format (see `dumpPoints` and `dumpMeshes` to know what format is)
     * @return formatted string
     */
    auto GridQuadLinear::dump() const -> std::string {
        auto oss = ostringstream();
        oss << dumpPoints() << "\n";
        oss << dumpMeshes();
        return oss.str();
    }

}
