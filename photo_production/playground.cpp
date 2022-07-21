#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

int main(){
    std::string const file = "C:/Users/deniz/Documents/GitHub/feynumeric/photo_production/data/experimental_data_pi_p_elastic.txt";

    std::ifstream in(file);

    std::string buffer;
    buffer.reserve(100UL);

    struct Row{
        double srt, cross, dcros, plab;
    };

    while( std::getline(in, buffer) ){
        if( buffer.starts_with('#') ){
            continue;
        }
        std::stringstream stream(buffer);
        Row row;
        stream >> row.srt >> row.cross >> row.dcros >> row.plab;
        std::cout << row.srt << " " << row.cross << " " << row.dcros << " " << row.plab << "\n";
    }
}