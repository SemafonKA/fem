#pragma once
#include <chrono>

class Timer {
public:
    inline Timer() {
        start();
    }

    void inline start() {
        startTime = std::chrono::system_clock::now();
        isRunning = true;
    }

    inline void stop() {
        endTime = std::chrono::system_clock::now();
        isRunning = false;
    }

    [[nodiscard]]
    inline auto elapsedMilliseconds() const -> int64_t {
        auto end = isRunning ? std::chrono::system_clock::now() : endTime;

        return std::chrono::duration_cast<std::chrono::milliseconds>(end - startTime).count();
    }

    [[nodiscard]]
    inline auto elapsedSeconds() const -> double {
        return static_cast<double>(elapsedMilliseconds()) / 1000.0;
    }

private:
    std::chrono::time_point<std::chrono::system_clock> startTime;
    std::chrono::time_point<std::chrono::system_clock> endTime;

    bool isRunning = false;
};