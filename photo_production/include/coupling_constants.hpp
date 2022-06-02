#ifndef FEYNUMERIC_COUPLING_CONSTANTS_HPP
#define FEYNUMERIC_COUPLING_CONSTANTS_HPP

#include <array>
#include <string>


inline std::array<std::string, 3> sort(std::string a, std::string b, std::string c){
    std::array<std::string, 3> arr{a,b,c};
    std::sort(arr.begin(), arr.end());
    return arr;
}

inline std::string coupling_string(std::string a, std::string b, std::string c){
    auto arr = sort(a, b, c);
    return arr[0]+arr[1]+arr[2];
}

#endif //FEYNUMERIC_COUPLING_CONSTANTS_HPP
