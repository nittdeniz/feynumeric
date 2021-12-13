#ifndef Feynumeric_TABLE_HPP
#define Feynumeric_TABLE_HPP

#include <functional>
#include <string>
#include <vector>

#include "range.hpp"

namespace Feynumeric
{
    using std::vector, std::string;
    void table(vector<string> const& table_header, std::function<vector<double>(double)> const& data, Range const& range);
}

#endif // Feynumeric_TABLE_HPP