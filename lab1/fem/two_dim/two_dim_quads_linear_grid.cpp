#include "two_dim_quads_linear_grid.h"

#include <vector>
#include <format>
#include <string>
#include <array>

#include "../../timer.h"
#include "../../logger.h"

using std::vector;
using std::array;
using std::format;
using std::string;
using fem::two_dim::Domain;

/**
 * @brief Returns Kx and Ky for grid (count of points by X and Y sides)
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
 */
static void fillSubdivides(fem::two_dim::GridQuadLinear& grid, size_t beg, size_t end, size_t step, double cx) {
    auto xLenght = grid.points.at(end).x - grid.points.at(beg).x;
    auto yLenght = grid.points.at(end).y - grid.points.at(beg).y;
    size_t nx = (end - beg) / step;
    auto fraction = cx == 1.0 ? (1.0 / static_cast<double>(nx)) : ((cx - 1.0) / (std::pow(cx, nx) - 1));
    auto xh1 = xLenght * fraction;
    auto yh1 = yLenght * fraction;

    for (size_t i = 1; i < nx; i++) {
        grid.points.at(beg + i * step).x = grid.points.at(beg).x + xh1 * std::pow(cx, i - 1);
        grid.points.at(beg + i * step).y = grid.points.at(beg).y + yh1 * std::pow(cx, i - 1);
    }
}

/**
 * @brief Compute and fill point coordinates for chosen domain
 */
static void fillPoints(const Domain& domain, fem::two_dim::GridQuadLinear& grid) {
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

namespace fem::two_dim {

    void GridQuadLinear::buildFrom(const Domain& domain) {
        // Get overall point count
        auto [countX, countY] = getPointsCount(domain);
        Kx = countX; Ky = countY;
        points.resize(Kx * Ky);
        logger::debug(format("Grid will be {}x{} size ({} overall)", Kx, Ky, points.size()));

        // Fill points
        auto timer = Timer();
        fillPoints(domain, *this);
        timer.stop();
        logger::debug(format("Points was computed in {} ms", timer.elapsedMilliseconds()));
    }

}
