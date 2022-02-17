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

	Feynman_Diagram_Ptr s_channel = std::make_shared<Feynman_Diagram>(Double_Wrench, VMP,
	                                                        std::initializer_list<Particle_Ptr>{Electron, Positron},
	                                                        std::initializer_list<Particle_Ptr>{Photon},
	                                                        std::initializer_list<Particle_Ptr>{Electron, Positron});

    Feynman_Process e_scattering({s_channel});
    e_scattering.dsigma_dcos_table(std::cout, 500._MeV, 0.1);
    return EXIT_SUCCESS;
}

