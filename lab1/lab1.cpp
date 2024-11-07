#include "logger.h"

using namespace std;

int main() {
    setlocale(LC_ALL, "ru-RU");
    logger::inFrameDebug("Включен режим отладки. Производительность будет снижена");

    return 0;
}
