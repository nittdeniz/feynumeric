#include "feynman_graph.hpp"
#include "matrix.hpp"
#include "particles.hpp"
#include "qed.hpp"
#include "vertex.hpp"

namespace Feynumeric
{
	Vertex electron_electron_photon({
			                                {Electron, Direction::INCOMING | Direction::OUTGOING},
			                                {Electron, Direction::INCOMING | Direction::OUTGOING},
			                                {Photon, Direction::INCOMING | Direction::OUTGOING},
	                                },
								 [](Kinematics const& kin, Particle_List const& particles)
								 {
									return Matrix(1,1,kin.sqrt_s() * particles[0]->spin()->j());
								 }
								 );
	Vertex electron_positron_photon({
			                                {Electron, Direction::INCOMING | Direction::OUTGOING},
			                                {Positron, Direction::INCOMING | Direction::OUTGOING},
			                                {Photon, Direction::INCOMING | Direction::OUTGOING},
	                                },
	                                [](Kinematics const& kin, Particle_List const& particles)
	                                {
		                                return Matrix(1,1,kin.sqrt_s() * particles[0]->spin()->j());
	                                });
}
