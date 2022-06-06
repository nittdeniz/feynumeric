#ifndef FEYNUMERIC_COUPLINGS_HPP
#define FEYNUMERIC_COUPLINGS_HPP

#include <map>
#include <string>

class Couplings
{
private:
    std::map<std::string, double> _constants;
public:
    Couplings(std::string const& file_in);
    double get(std::string const& key) const;
    void set(std::string const& key, double value);
};
inline std::array<std::string, 3> sort(std::string a, std::string b, std::string c){
    std::array<std::string, 3> arr{a,b,c};
    std::sort(arr.begin(), arr.end());
    return arr;
}

inline std::string coupling_string(std::string a, std::string b, std::string c){
    auto arr = sort(a, b, c);
    return "g"+arr[0]+arr[1]+arr[2];
}


#endif //FEYNUMERIC_COUPLINGS_HPP
