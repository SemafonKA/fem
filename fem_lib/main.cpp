#include <fstream>
#include <iostream>
#include <sstream>
#include <numbers>

#include "logger.h"
#include "fem/includes/fem.h"


using namespace std;
using namespace fem::three_dim;


//double u(double x, double y, double z) {
//    return x + y + z;
//}

double lambda([[maybe_unused]] double x, [[maybe_unused]] double y, [[maybe_unused]] double z, [[maybe_unused]] size_t material) {
    return 1;
}

double gamma([[maybe_unused]] double x, [[maybe_unused]] double y, [[maybe_unused]] double z, [[maybe_unused]] size_t material) {
    return 1;
}

//double func([[maybe_unused]] double x, [[maybe_unused]] double y, [[maybe_unused]] double z, [[maybe_unused]] size_t material) {
//    return gamma(x, y, z, material) * u(x, y, z);
//}
//
//double s1_func([[maybe_unused]] double x, [[maybe_unused]] double y, [[maybe_unused]] double z, [[maybe_unused]] size_t index) {
//    return u(x, y, z);
//}
//
//double s2_func([[maybe_unused]] double x, [[maybe_unused]] double y, [[maybe_unused]] double z, [[maybe_unused]] size_t index) {
//    return DomainFunctionsCuboidLinear::not_defined;
//}


double u_t(double x, double y, double z, double t) {
    return sin(t);
}

double func_t(double x, double y, double z, double t, size_t material) {
    return cos(t);
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

    auto funcs = DomainFunctionsCuboidLinearDynamic();
    funcs.sigma = gamma;
    funcs.lambda = lambda;
    funcs.func = func_t;
    funcs.s1_func = s1_func_t;
    funcs.initials = u_t;

    double t_beg = 0.0;
    double t_end = 2.0;
    size_t t_count = 3;
    double t_step = (t_end - t_beg) / (t_count - 1);
    vector<double> meshT(t_count);
    meshT[0] = t_beg;
    meshT[t_count - 1] = t_end;
    for (size_t i = 1; i < t_count - 1; i++) {
        meshT[i] = t_beg + t_step * i;
    }

    auto solver = SolverCuboidLinear(funcs, grid);
    auto q = solver.solveStatic(meshT);

    return 0;
}
