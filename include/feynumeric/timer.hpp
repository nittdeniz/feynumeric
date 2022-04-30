#ifndef FEYNUMERIC_TIMER_HPP
#define FEYNUMERIC_TIMER_HPP
#include <chrono>

namespace Feynumeric
{
    class Timer{
    private:
        bool _is_running;
        std::chrono::time_point <std::chrono::steady_clock> _start_time;
        std::chrono::time_point <std::chrono::steady_clock> _end_time;
    public:
        void start();
        void stop();
        template <typename T> double time()
        {
            std::chrono::time_point <std::chrono::steady_clock> endTime;
            auto end_time = _is_running? std::chrono::steady_clock::now() : _end_time;
            return std::chrono::duration_cast<T>(end_time - _start_time).count();
        }
    };
}
#endif