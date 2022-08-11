#include "format.hpp"
#include "messages.hpp"
#include "table.hpp"

#include <fstream>
#include <sstream>
#include <string>

namespace Feynumeric{
	Table::Table(std::map<double, double> const& values)
	: _values(values)
    , lower_out_of_bounds(false)
    , upper_out_of_bounds(false)
	{

	}

	double Table::interpolate(double value) const{
		if( _values.contains(value) ){
			return _values.at(value);
		}
		auto upper_bound = _values.upper_bound(value);
//        auto lower_bound = _values.lower_bound(value);
//		if( !upper_out_of_bounds && upper_bound == _values.end() ){
//			warning(FORMAT("Value {} is out of range. Can not interpolate and will return last value.", value));
//            upper_out_of_bounds = true;
//		}
//		if(  !lower_out_of_bounds && lower_bound == _values.begin() ){
//			warning(FORMAT("Value {} is out of range. Can not interpolate and will return first value.", value));
//            lower_out_of_bounds = true;
//		}
		auto const a = --upper_bound;
		auto const b = ++upper_bound;
		return a->second + (b->second-a->second)/(b->first-a->first) * (value - a->first);
	}

	std::ofstream& operator<<(std::ofstream& out, Table const& table){
		for( auto const& [key, value] : table._values ){
			out << FORMAT("{} {}\n", key, value);
		}
		return out;
	}

	std::ifstream& operator>>(std::ifstream& in, Table& table){
		std::string buffer;
		while( std::getline(in, buffer) ){
			if( buffer.starts_with('#') ){
				continue;
			}
			std::stringstream stream(buffer);
			double key, value;
			stream >> key >> value;
			table._values[key] = value;
		}
		return in;
	}

	void Table::write(std::string file_out){
		std::ofstream out(file_out);
		if( out ){
			out << *this;
		}
	}

	void Table::read(std::string file_in){
		std::ifstream in(file_in);
		if( in ){
			in >> *this;
		}
	}
}

