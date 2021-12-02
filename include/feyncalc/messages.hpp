#ifndef FEYNCALC_MESSAGES_HPP
#define FEYNCALC_MESSAGES_HPP

#include <string>

namespace Feyncalc
{
    [[noreturn]] void critical_error(std::string const& message);
    void error(std::string const& message);
    void warning(std::string const& message);
    void status(std::string const& message);
}



#endif // FEYNCALC_MESSAGES_HPP