#include "write_map.hpp"

#include <fstream>
#include <string>

void write_map(const std::map<double, double> &map, std::string const& file_name){
    std::ofstream out(file_name);
    if( out ){
        for( auto const&[key, value] : map )
        {
            out << key << ", " << value << "\n";
        }
    }
    out.close();
}
