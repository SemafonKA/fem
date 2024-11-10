#include <fstream>

#include "logger.h"
#include "fem/includes/fem.h"

using namespace std;
using namespace fem::two_dim;

int main() {
    setlocale(LC_ALL, "ru-RU");

    const auto domainFilepath = string("domain.txt");
    const auto gridFilepath = string("grid.txt");

    logger::inFrameDebug("Debug mode enabled. Program may work slow and additional logs was output");

    logger::log(format("Reading domain from file \"{}\"", domainFilepath));
    Domain domain;
    try {
        // TODO: ѕроработать обработку ошибок при чтении домена
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

    return 0;
}
