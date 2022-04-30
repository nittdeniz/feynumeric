#include "timer.hpp"

namespace Feynumeric{
    void Timer::start(){
        _start_time = std::chrono::steady_clock::now();
        _is_running = true;
    }

    void Timer::stop()
    {
        _end_time = std::chrono::steady_clock::now();
        _is_running = false;
    }
}
