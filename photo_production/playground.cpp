#include <feynumeric/matrix.hpp>
#include <feynumeric/polynomial.hpp>

#include <iostream>

int main(){
	using namespace Feynumeric;

	std::vector<std::pair<double, Complex>> list = {{1., 3.75872}, {2., 9.46583}, {3., 25.8396}, {4., 58.88}, {5.,
	                                                    114.587}};

	Polynomial f(1);
	f.fit(list);

}