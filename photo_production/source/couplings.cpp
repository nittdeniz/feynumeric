#include "couplings.hpp"

#include <feynumeric/format.hpp>
#include <feynumeric/messages.hpp>

#include <fstream>
#include <sstream>

Couplings::Couplings(const std::string &file_in)
{
    std::ifstream ifs(file_in);
    if( !ifs ){
        Feynumeric::critical_error(FORMAT("Could not open {}.", file_in));
    }

    std::string buffer;
    while( std::getline(ifs, buffer) )
    {
        std::string key;
        double value;
        if( buffer.starts_with('#'))
        {
            continue;
        }

        std::stringstream sstream(buffer);
        sstream >> key >> value;
        _constants[key] = value;
    }
}

double Couplings::get(const std::string &key) const
{
    if( !_constants.contains(key) ){
        Feynumeric::warning(FORMAT("Coupling constant {} does not exist. Returning 1.\n", key));
        return 1.;
    }
    return _constants.at(key);
}

void Couplings::set(const std::string &key, double value)
{
    _constants[key] = value;
}
