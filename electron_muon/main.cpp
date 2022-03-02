#include <initializer_list>
#include <iostream>
#include <memory>

#include <feynumeric/feynman_diagram.hpp>
#include <feynumeric/feynman_process.hpp>
#include <feynumeric/units.hpp>

#include <feynumeric/dirac.hpp>
#include <feynumeric/qed.hpp>

#include <feynumeric/topologies.hpp>


int main()
{
    using namespace Feynumeric;
    using namespace Feynumeric::Units;
    using namespace QED;

    init_particles();
	init_vertices();

	auto s_channel = create_diagram("s_channel", Double_Wrench, VMP,
	                                {Electron, Positron},
	                                {Photon},
	                                {Electron, Positron});

	double const q = momentum(0.5, Electron->mass(), Electron->mass());

	double const cos = 0.9125004096970463;

	auto s1 = std::make_shared<Angular_Momentum>(0.5, 0.5);

	std::cerr << ubar(Electron, four_momentum(q, Electron->mass(), cos), s1, {}) << "\n";
	std::cerr << v(Positron, four_momentum(-q, Positron->mass(), cos), s1, {}) << "\n";

	Feynman_Process pair_production({s_channel});
	pair_production.dsigma_dcos_table(std::cout, 500._MeV, {{cos}});

//	Feynman_Diagram_Ptr s_channel = create_diagram("s_channel", Double_Wrench, VMP,
//	                                               {Electron, Photon},
//	                                               {Electron},
//	                                               {Photon, Electron});
//
//	Feynman_Diagram_Ptr t_channel = create_diagram("t_channel", X_Man, VMP,
//	                                               {Electron, Photon},
//	                                               {Electron},
//	                                               {Photon, Electron});


//	Feynman_Process compton({s_channel, t_channel});

//	std::stringstream out;

//	double const cos_theta = 0.2134;

//	compton.dsigma_dcos_table(std::cout, 500._MeV, 0.1);
    return EXIT_SUCCESS;
}

