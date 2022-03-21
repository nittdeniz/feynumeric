#include <feynumeric/feynumeric.hpp>

#include <iostream>

int main(int argc, char** argv)
{
	std::string file_name(argv[1]);
	volatile Feynumeric::Particle_Manager PM(file_name);
	std::cout << "Hello World\n";
}

