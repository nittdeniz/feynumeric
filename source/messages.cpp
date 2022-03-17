#include "messages.hpp"

#include <iostream>

namespace Feynumeric
{
    [[noreturn]] void critical_error(const std::string &message)
    {
        std::cerr << "!!! Critical Error: " << message << "\n";
        std::abort();
    }

    void error(const std::string &message)
    {
        std::cerr << "@ Error: " << message << "\n";
    }

    void warning(const std::string &message)
    {
        std::cerr << "$ Warning: " << message << "\n";
    }

    void status(const std::string &message)
    {
        std::cerr << "[Status] " << message << "\n";
    }

	critical_exception::critical_exception(std::string&& message)
	: _message(std::move(message))
	{

	}
}

