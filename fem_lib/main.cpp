#include <fstream>
#include <iostream>
#include <sstream>
#include <numbers>

#include "logger.h"
#include "fem/includes/fem.h"


using namespace std;
using namespace fem::three_dim;


double u(double x, double y, double z) {
    return x + y + z;
}

double lambda([[maybe_unused]] double x, [[maybe_unused]] double y, [[maybe_unused]] double z, [[maybe_unused]] size_t material) {
    return 1;
}

double gamma([[maybe_unused]] double x, [[maybe_unused]] double y, [[maybe_unused]] double z, [[maybe_unused]] size_t material) {
    return 1;
}

double func([[maybe_unused]] double x, [[maybe_unused]] double y, [[maybe_unused]] double z, [[maybe_unused]] size_t material) {
    return gamma(x, y, z, material) * u(x, y, z);
}

double s1_func([[maybe_unused]] double x, [[maybe_unused]] double y, [[maybe_unused]] double z, [[maybe_unused]] size_t index) {
    return u(x, y, z);
}

double s2_func([[maybe_unused]] double x, [[maybe_unused]] double y, [[maybe_unused]] double z, [[maybe_unused]] size_t index) {
    return DomainFunctionsCuboidLinear::not_defined;
}


double u_t(double x, double y, double z, double t) {
    return x + y + z;
}

double func_t(double x, double y, double z, double t, size_t material) {
    return gamma(x, y, z, material) * u_t(x, y, z, t);
}

double s1_func_t(double x, double y, double z, double t, size_t index) {
    return u_t(x, y, z, t);
}



void printSolveToFile(const string& filepath, size_t pxCount, size_t pyCount, const SolverCuboidLinear& solver) {
    auto ofile = ofstream();
    ofile.open(filepath);

    auto& grid = solver.getGrid();
    if (pxCount == 0 && pyCount == 0) {
        auto& solution = solver.getSolution();
        for (size_t i = 0; i < grid.points.size(); i++) {
            auto& p = grid.points[i];
            ofile << format("{:20.14f} {:20.14f} {:20.14f}\n", p.x, p.y, solution.at(i));
        }
    }
    else {
        // ??? TODO: Добавить вывод произвольной области, ы
    }

    ofile.close();
}

int main() {
    setlocale(LC_ALL, "ru-RU");

    const auto domainFilepath = string("io/domain.txt");
    const auto gridFilepath = string("io/grid.txt");
    const auto solveFilepath = string("io/solution.txt");

    logger::inFrameDebug("Debug mode enabled. Program may work slow and additional logs was output");

    logger::log(format("Reading domain from file \"{}\"", domainFilepath));
    Domain domain;
    try {
        // TODO: Проработать обработку ошибок при чтении домена
        domain = Domain::readFromFile(domainFilepath);
    }
    catch (const std::runtime_error& e) {
        logger::error(e.what());
        logger::error("This error break all further program logic so the program will be terminated");
        return -1;
    }
    logger::log("Reading domain from file was successful");

    auto grid = GridCuboidLinear();
    try {
        grid.buildFrom(domain);
    }
    catch (const std::runtime_error& e) {
        logger::error(format("There is some error while building grid: {}", e.what()));
        return -1;
    }

    logger::log("Grid was builded");
    logger::log(format("Count of nodes: {}, count of elements: {}", grid.points.size(), grid.meshes.size()));

    logger::debug("Finite elements dump:");
    for (const auto& elem : grid.meshes) {
        logger::debug(format("   {}", elem.toString()));
    }

    auto ofile = ofstream();
    ofile.open(gridFilepath);
    ofile << grid.dump();
    ofile.close();
    logger::log(format("Grid info was written into file {}", gridFilepath));

    auto funcs = DomainFunctionsCuboidLinear();
    funcs.gamma = gamma;
    funcs.lambda = lambda;
    funcs.func = func;
    funcs.s1_func = s1_func;
    funcs.s2_func = s2_func;

    auto solver = SolverCuboidLinear(funcs, grid);
    auto q = solver.solveStatic();

    auto absolute = 0.0;
    auto relative = 0.0;

    logger::inFrame("Solver result");
    auto oss = ostringstream();

    if (grid.points.size() < 100) {
        oss << format("| {:^4} | {:^14} | {:^14} | {:^14} | {:^14} | {:^14} | {:^14} | {:^14} |\n",
            "#", "X", "Y", "Z", "q_i", "u_i", "abs", "rel");
        oss << format("| {0:^4} | {0:^14} | {0:^14} | {0:^14} | {0:^14} | {0:^14} | {0:^14} | {0:^14} |\n", ":-:");
    }

    for (size_t i = 0; i < grid.points.size(); i++) {
        const auto& point = grid.points[i];
        auto qi = q[i];
        auto ui = u(point.x, point.y, point.z);
        auto p_abs = std::abs(qi - ui);
        auto p_rel = p_abs / std::abs(ui);

        absolute += p_abs;
        relative += std::abs(ui); // Накапливаем значение общей функции

        if (grid.points.size() < 100) {
            oss << format("| {:^4} | {:^14f} | {:^14f} | {:^14f} | {:^14f} | {:^14f} | {:^14f} | {:^14f} |\n",
                i + 1, point.x, point.y, point.z, qi, ui, p_abs, p_rel);
        }
    }

    cout << oss.str() << "\n";
    relative = absolute / relative;
    cout << format("Absolute error: {:8.4e}\n", absolute);
    cout << format("Relative error: {:8.4e}\n", relative);

    logger::inFrame("Test points");
    std::vector<Point> testPoints = {
        {2.5, 1.4, 2.2},
    };

    oss.str("");
    relative = absolute = 0;

    oss << format("| {:^4} | {:^14} | {:^14} | {:^14} | {:^14} | {:^14} | {:^14} | {:^14} |\n",
        "#", "X", "Y", "Z", "q_i", "u_i", "abs", "rel");
    oss << format("| {0:^4} | {0:^14} | {0:^14} | {0:^14} | {0:^14} | {0:^14} | {0:^14} | {0:^14} |\n", ":-:");

    for (size_t i = 0; i < testPoints.size(); i++) {
        const auto& point = testPoints[i];
        auto qi = solver.value(point);
        auto ui = u(point.x, point.y, point.z);
        auto p_abs = std::abs(qi - ui);
        auto p_rel = p_abs / std::abs(ui);

        absolute += p_abs;
        relative += std::abs(ui); // Накапливаем значение общей функции

        oss << format("| {:^4} | {:^14f} | {:^14f} | {:^14f} | {:^14f} | {:^14f} | {:^14f} | {:^14f} |\n",
            i + 1, point.x, point.y, point.z, qi, ui, p_abs, p_rel);
    }
    cout << oss.str() << "\n";
    relative = absolute / relative;
    cout << format("Absolute error: {:8.4e}\n", absolute);
    cout << format("Relative error: {:8.4e}\n", relative);

    printSolveToFile(solveFilepath, 0, 0, solver);

    return 0;
}
