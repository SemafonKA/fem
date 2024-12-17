#include <fstream>
#include <iostream>
#include <sstream>

#include "logger.h"
#include "fem/includes/fem.h"


using namespace std;
using namespace fem::two_dim;


double u(double x, double y) {
    return x + y;
}

double lambda([[maybe_unused]] double x, [[maybe_unused]] double y, [[maybe_unused]] size_t material) {
    return 1.0;
}

double gamma([[maybe_unused]] double x, [[maybe_unused]] double y, [[maybe_unused]] size_t material) {
    return 0.01;
}

double func([[maybe_unused]] double x, [[maybe_unused]] double y, [[maybe_unused]] size_t material) {
    return gamma(x, y, material) * u(x, y);
}

double s1_func([[maybe_unused]] double x, [[maybe_unused]] double y, [[maybe_unused]] size_t index) {
    return u(x, y);
}

double s2_func([[maybe_unused]] double x, [[maybe_unused]] double y, [[maybe_unused]] size_t index) {
    return DomainFunctionsQuadsLinear::not_defined;
}

int main() {
    setlocale(LC_ALL, "ru-RU");

    const auto domainFilepath = string("io/domain.txt");
    const auto gridFilepath = string("io/grid.txt");

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

    auto grid = GridQuadLinear();
    try {
        grid.buildFrom(domain);
    }
    catch (const std::runtime_error& e) {
        logger::error(format("There is some error while building grid: {}", e.what()));
        return -1;
    }
    auto ofile = ofstream();
    ofile.open(gridFilepath);
    ofile << grid.dump();
    ofile.close();
    logger::log(format("Grid info was written into file {}", gridFilepath));
    
    auto funcs = DomainFunctionsQuadsLinear();
    funcs.gamma = gamma;
    funcs.lambda = lambda;
    funcs.func = func;
    funcs.s1_func = s1_func;
    funcs.s2_func = s2_func;

    auto solver = SolverQuadsLinear(funcs, grid);
    auto q = solver.solveStatic();

    auto absolute = 0.0;
    auto relative = 0.0;

    logger::inFrame("Solver result");
    auto oss = ostringstream();

    if (grid.points.size() < 100) {
        oss << format("| {:^4} | {:^22} | {:^22} | {:^22} | {:^22} | {:^22} | {:^22} |\n",
            "#", "X", "Y", "q_i", "u_i", "abs", "rel");
        oss << format("| {0:^4} | {0:^22} | {0:^22} | {0:^22} | {0:^22} | {0:^22} | {0:^22} |\n", ":-:");
    }

    for (size_t i = 0; i < grid.points.size(); i++) {
        const auto& point = grid.points[i];
        auto qi = q[i];
        auto ui = u(point.x, point.y);
        auto p_abs = std::abs(qi - ui);
        auto p_rel = p_abs / std::abs(ui);

        absolute += p_abs;
        relative += std::abs(ui); // Накапливаем погрешность

        if (grid.points.size() < 100) {
            oss << format("| {:^4} | {:^22f} | {:^22f} | {:^22f} | {:^22f} | {:^22f} | {:^22f} |\n",
                i + 1, point.x, point.y, qi, ui, p_abs, p_rel);
        }
    }

    cout << oss.str() << "\n";
    relative = absolute / relative;
    cout << format("Absolute error: {:8.4e}\n", absolute);
    cout << format("Relative error: {:8.4e}\n", relative);

    logger::inFrame("Test points");
    std::vector<Point> testPoints = {
        {1, 1},
        {1.25, 1.5},
        {5.2, 1.444},
        {3.2, 2.5},
        {2.8, 3.4}
    };

    oss.str("");
    relative = absolute = 0;

    oss << format("| {:^4} | {:^22} | {:^22} | {:^22} | {:^22} | {:^22} | {:^22} |\n",
        "#", "X", "Y", "q_i", "u_i", "abs", "rel");
    oss << format("| {0:^4} | {0:^22} | {0:^22} | {0:^22} | {0:^22} | {0:^22} | {0:^22} |\n", ":-:");

    for (size_t i = 0; i < testPoints.size(); i++) {
        const auto& point = testPoints[i];
        auto qi = solver.value(point);
        auto ui = u(point.x, point.y);
        auto p_abs = std::abs(qi - ui);
        auto p_rel = p_abs / std::abs(ui);

        absolute += p_abs;
        relative += std::abs(ui); // Накапливаем погрешность

        oss << format("| {:^4} | {:^22f} | {:^22f} | {:^22f} | {:^22f} | {:^22f} | {:^22f} |\n",
            i + 1, point.x, point.y, qi, ui, p_abs, p_rel);
    }
    cout << oss.str() << "\n";
    relative = absolute / relative;
    cout << format("Absolute error: {:8.4e}\n", absolute);
    cout << format("Relative error: {:8.4e}\n", relative);

    return 0;
}
