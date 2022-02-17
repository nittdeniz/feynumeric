#include <initializer_list>
#include <iostream>
#include <memory>

#include <feynumeric/feynman_diagram.hpp>
#include <feynumeric/feynman_process.hpp>
#include <feynumeric/process.hpp>
#include <feynumeric/topology.hpp>
#include <feynumeric/units.hpp>
#include <feynumeric/qed.hpp>


#include "particles.hpp"
#include "vertices.hpp"

#include <feynumeric/dirac.hpp>
#include <feynumeric/topologies.hpp>


int main()
{
    using namespace Feynumeric;
    using namespace Feynumeric::Units;

    init_particles();
	init_vertices();

//	Feynman_Diagram_Ptr s_channel = std::make_shared<Feynman_Diagram>(
//								Double_Wrench,
//								VMP,
//	                          std::initializer_list<Particle_Ptr>({Photon, Electron}),
//	                          std::initializer_list<Particle_Ptr>({Electron}),
//	                          std::initializer_list<Particle_Ptr>({Electron, Photon})
//	                          );

	Feynman_Diagram_Ptr t_channel = std::make_shared<Feynman_Diagram>(X_Man, VMP,
	                          std::initializer_list<Particle_Ptr>({Photon, Electron}),
	                          std::initializer_list<Particle_Ptr>({Positron}),
	                          std::initializer_list<Particle_Ptr>({Electron, Photon})
	);

    Feynman_Process e_muon_scattering({t_channel});
    e_muon_scattering.dsigma_dcos_table(std::cout, 2._MeV, 20ULL);
//    return EXIT_SUCCESS;
}

