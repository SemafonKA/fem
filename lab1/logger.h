#pragma once
#include <iostream>
#include <format>
#include <chrono>
#include <sstream>


#include "termcolor/termcolor.h"

namespace logger {

    enum class Colors {
        red,
        green,
        debug_green,
        blue,
        yellow,
        standard,
        success,
        error,
    };

    /**
     * @brief Logs simple text into terminal with timestamp added
     * @param str - text to log into terminal (without `\n` symbol)
     * @param color - color to be outputted to the terminal
     */
    void log(const std::string& str, Colors color = Colors::standard);

    /**
     * @brief Logs error text into terminal with timestamp added
     * @param str - text of error to log into terminal (without `\n` symbol)
     * @param color - color to be outputted to the terminal
     */
    void error(const std::string& str, Colors color = Colors::red);

    /**
     * @brief Logs simple debug text into terminal with timestamp added. Work ONLY in debug mode!
     * @param str - text to log into terminal (without `\n` symbol)
     * @param color - color to be outputted to the terminal (debug_green by default)
     */
    void debug(const std::string& str, Colors color = Colors::debug_green);

    /**
     * @brief Logs some text into terminal into starred (`*`) frame
     * @param str - text to log into frame (should be not so long!)
     * @param color - color to be outputted to the terminal
     */
    void inFrame(const std::string& str, Colors color = Colors::standard);

    /**
     * @brief Logs some debug text into terminal into starred (`*`) frame
     * @param str - text to log into frame (should be not so long!)
     * @param color - color to be outputted to the terminal (debug_green by default)
     */
    void inFrameDebug(const std::string& str, Colors color = Colors::debug_green);
}
