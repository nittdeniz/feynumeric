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

	auto s_channel = create_diagram("t_channel", X_Man, VMP,
	                                {Photon, Electron},
	                                {Positron},
	                                {Electron, Photon});


	Feynman_Process pair_production({s_channel});
	pair_production.dsigma_dcos_table(std::cout, 500._MeV, 0.1);
	compton.dsigma_dcos_table(std::cout, 500._MeV, 0.1);
    return EXIT_SUCCESS;
}

