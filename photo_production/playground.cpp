#include <feynumeric/feynumeric.hpp>
#include <feynumeric/qed.hpp>

#include <iostream>

int main(int argc, char** argv)
{
	using namespace Feynumeric;
	std::string file_name(argv[1]);
	volatile Particle_Manager PM(file_name);

	Vertex_Manager VM;

	VM.add(Vertex(
			{
					{QED::Positron, Direction::ANY},
					{QED::Electron, Direction::ANY},
					{QED::Photon, Direction::OUTGOING}
		},
		[](Kinematics const&, Particle_List const& list){
				std::cout << list[0]->particle()->name() << "\n";
				std::cout << list[1]->particle()->name() << "\n";
				std::cout << list[2]->particle()->name() << "\n";
				return Matrix(1,1,1);}
	));

	VM.add(Vertex(
			{
					{QED::Electron, Direction::ANY},
					{QED::Electron, Direction::ANY},
					{QED::Photon, Direction::OUTGOING}
			},
			[](Kinematics const&, Particle_List const&){return Matrix(1,1,1);}
	));

	VM.find({QED::Photon, Direction::OUTGOING}
			{QED::Positron, Direction::INCOMING},
	        {QED::Electron, Direction::INCOMING},
	        );

	std::cout << "end\n";

}

