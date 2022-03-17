#ifndef Feynumeric_MESSAGES_HPP
#define Feynumeric_MESSAGES_HPP

#include <string>

#ifdef TEST_MODE
	#define CRITICAL_ERROR(msg) throw critical_exception(msg);
#else
	#define CRITICAL_ERROR(msg) critical_error(msg);
#endif

namespace Feynumeric
{
    [[noreturn]] void critical_error(std::string const& message);
    void error(std::string const& message);
    void warning(std::string const& message);
    void status(std::string const& message);

    class critical_exception : public std::exception{
    	std::string _message;
    public:
    	critical_exception(std::string&& message);
    };
}



#endif // Feynumeric_MESSAGES_HPP