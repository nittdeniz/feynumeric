#include <initializer_list>
#include <iostream>
#include <memory>

#include <feynumeric/feynman_diagram.hpp>
#include <feynumeric/feynman_process.hpp>
#include <feynumeric/units.hpp>

#include <feynumeric/qed.hpp>

#include <feynumeric/topologies.hpp>


int main()
{
    using namespace Feynumeric;
    using namespace Feynumeric::Units;
    using namespace QED;

    init_particles();
	init_vertices();

	Feynman_Diagram_Ptr s_channel = create_diagram("s_channel", Double_Wrench, VMP,
	                                               {Electron, Positron},
	                                               {Photon},
	                                               {Electron, Positron});

	Feynman_Diagram_Ptr t_channel = create_diagram("t_channel", X_Man, VMP,
	                                               {Electron, Positron},
	                                               {Photon},
	                                               {Electron, Positron});


	Feynman_Process e_scattering({s_channel, t_channel});

	std::stringstream out;

//	double const cos_theta = 0.2134;

	auto x =e_scattering.dsigma_dcos_table(500._MeV, {{0.1}});
    return EXIT_SUCCESS;
}

