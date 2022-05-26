#ifndef Feynumeric_TABLE_HPP
#define Feynumeric_TABLE_HPP

#include <fstream>
#include <map>

namespace Feynumeric{
	class Table{
		std::map<double, double> _values;
	public:
		Table(std::map<double, double> const& values);
		double interpolate(double value) const;
		friend std::ifstream& operator>>(std::ifstream& in, Table& table);
		friend std::ofstream& operator<<(std::ofstream& out, Table const& table);

		void write(std::string file_out);
		void read(std::string file_in);
	};
}

#endif // Feynumeric_TABLE_HPP