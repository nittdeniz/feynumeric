#ifndef Feynumeric_TABLE_HPP
#define Feynumeric_TABLE_HPP
/*
#include <functional>
#include <string>
#include <vector>
#include "format.hpp"


namespace Feynumeric
{
	template<typename T>
    void table(std::ostream& out, std::vector<std::string> const& table_header, std::vector<std::function<std::vector<T>(int)>> const& data, int const start, int const end, int const step_size = 1)
	{
    	for( auto const& header : table_header )
	    {
		    out << FORMAT("|{:^10}", header);
	    }
    	out << "|\n";
    	for( int a = start, i = 0; a <= end && i < data.size(); a += step_size, i++ )
	    {
    		for( auto const& item : data )
		    {
    			out << FORMAT("|{:^10}",
		    }
	    }
	}
}
*/
#endif // Feynumeric_TABLE_HPP