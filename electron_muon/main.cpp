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

	Feynman_Diagram_Ptr t_channel = create_diagram("t_channel",
			X_Man,
			VMP,
			{Electron, Electron},
			{Photon},
			{Electron, Electron});

	Feynman_Diagram_Ptr u_channel = create_diagram("u_channel",
			X_Man,
			VMP,
			{Electron, Electron},
			{Photon},
			{Electron, Electron});
	u_channel->cross_outgoing(0, 1);

    Feynman_Process e_scattering({t_channel, u_channel});
    e_scattering.dsigma_dcos_table(std::cout, 500._MeV, 0.1);
    return EXIT_SUCCESS;
}

