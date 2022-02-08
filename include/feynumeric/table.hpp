#ifndef Feynumeric_TABLE_HPP
#define Feynumeric_TABLE_HPP

#include <functional>
#include <string>
#include <vector>

#include "range.hpp"

namespace Feynumeric
{
    void table(std::vector<std::string> const& table_header, std::function<std::vector<double>(double)> const& data, Range const& range);
}

#endif // Feynumeric_TABLE_HPP