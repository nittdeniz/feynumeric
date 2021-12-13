#ifndef Feynumeric_MESSAGES_HPP
#define Feynumeric_MESSAGES_HPP

#include <string>

namespace Feynumeric
{
    [[noreturn]] void critical_error(std::string const& message);
    void error(std::string const& message);
    void warning(std::string const& message);
    void status(std::string const& message);
}



#endif // Feynumeric_MESSAGES_HPP