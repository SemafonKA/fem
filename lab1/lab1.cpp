#include "logger.h"

#include "fem/includes/domain.h"

using namespace std;
using namespace fem::two_dim;

int main() {
    setlocale(LC_ALL, "ru-RU");
    logger::inFrameDebug("Debug mode enabled. Program may work slow and additional logs was output");

    string domainFilepath = "grid.txt";

    logger::log(format("Reading domain from file \"{}\"", domainFilepath));
    
    auto domain = Domain();
    try {
        domain = Domain::readFromFile(domainFilepath);
    }
    catch (const std::runtime_error& _) {
        logger::error("This error break all further program logic so the program will be terminated");
        return -1;
    }

    logger::log("Reading domain from file was successed", logger::Colors::success);

    cout << domain.dump() << endl;
    cout << domain.dumpNoComments() << endl;

    return 0;
}
