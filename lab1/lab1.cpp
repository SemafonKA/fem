#include "logger.h"

#include "fem/two_dim/domain_reader.h"

using namespace std;
using namespace fem::two_dim;

int main() {
    setlocale(LC_ALL, "ru-RU");
    logger::inFrameDebug("Включен режим отладки. Производительность будет снижена");

    string domainFilepath = "grid.txt";

    logger::log(format("Reading domain from file \"{}\"", domainFilepath));
    auto reader = readDomainFromFile(domainFilepath);

    return 0;
}
