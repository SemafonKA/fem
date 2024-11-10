#include "logger.h"

namespace logger {

    void inline setColor(Colors color) {
        switch (color) {
        case Colors::red:
        case Colors::error:
            std::cout << termcolor::red;
            std::cerr << termcolor::red;
            return;

        case Colors::green:
        case Colors::success:
            std::cout << termcolor::bright_green;
            std::cerr << termcolor::bright_green;
            return;

        case Colors::debug_green:
            std::cout << termcolor::green;
            std::cerr << termcolor::green;
            return;

        case Colors::blue:
            std::cout << termcolor::bright_blue;
            std::cerr << termcolor::bright_blue;
            return;

        case Colors::warning:
        case Colors::yellow:
            std::cout << termcolor::yellow;
            std::cerr << termcolor::yellow;
            return;

        case Colors::standard:
            std::cout << termcolor::reset;
            std::cerr << termcolor::reset;
            return;

        default:
            throw std::logic_error("Used param is not implemented right now");
        }
    }

    static inline void resetColor() {
        setColor(Colors::standard);
    }

    static inline auto getCurrentTime() -> std::string {
        auto const time = std::chrono::current_zone()->to_local(std::chrono::system_clock::now());
        return std::format("{:%X}", time); // add %Y-%m-%d to print date also
    }



    /**
     * @brief Logs simple text into terminal with timestamp added
     * @param str - text to log into terminal (without `\n` symbol)
     * @param color - color to be outputted to the terminal
     */
    void log(const std::string& str, Colors color) {
        auto time = getCurrentTime();
        setColor(color);
        std::cout << std::format("[MSG : {}] - {}\n", time, str);
        resetColor();
    }

    /**
     * @brief Logs warning text into terminal with timestamp added
     * @param str - text to log into terminal (without `\n` symbol)
     * @param color - color to be outputted to the terminal
     */
    void warn(const std::string& str, Colors color) {
        auto time = getCurrentTime();
        setColor(color);
        std::cerr << std::format("[WRN : {}] - {}\n", time, str);
        resetColor();
    }


    /**
     * @brief Logs error text into terminal with timestamp added
     * @param str - text of error to log into terminal (without `\n` symbol)
     * @param color - color to be outputted to the terminal
     */
    void error(const std::string& str, Colors color) {
        auto time = getCurrentTime();
        setColor(color);
        std::cerr << std::format("[ERR : {}] - {}\n", time, str);
        resetColor();
    }

    /**
     * @brief Logs simple debug text into terminal with timestamp added
     * @param str - text to log into terminal (without `\n` symbol)
     * @param color - color to be outputted to the terminal (green by default)
     */
    void debug(const std::string& str, Colors color) {
#ifndef NDEBUG
        auto time = getCurrentTime();
        setColor(color);
        std::cout << std::format("[DBG : {}] - {}\n", time, str);
        resetColor();
#endif
    }

    /**
     * @brief Logs some text into terminal into starred (`*`) frame
     * @param str - text to log into frame (should be not so long!)
     * @param color - color to be outputted to the terminal
     */
    void inFrame(const std::string& str, Colors color) {
        auto len = str.length();
        auto oss = std::ostringstream();
        oss << std::format("\n{0:*^{1}}\n", "", len + (3 * 2));
        oss << std::format("** {} **\n", str);
        oss << std::format("{0:*^{1}}\n\n", "", len + (3 * 2));

        setColor(color);
        std::cout << oss.str();
        resetColor();
    }

    /**
     * @brief Logs some debug text into terminal into starred (`*`) frame
     * @param str - text to log into frame (should be not so long!)
     * @param color - color to be outputted to the terminal (green by default)
     */
    void inFrameDebug(const std::string& str, Colors color) {
#ifndef NDEBUG
        inFrame(str, color);
#endif
    }
}

