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

//	Feynman_Diagram_Ptr s_channel_1 = create_diagram("s_channel 1", Double_Wrench, VMP,
//	                                               {Photon, Electron},
//	                                               {Electron},
//	                                               {Electron, Photon});
//
//	Feynman_Diagram_Ptr s_channel_2 = create_diagram("s_channel 2", Double_Wrench, VMP,
//	                                               {Photon, Electron},
//	                                               {Electron},
//	                                               {Photon,Electron});
//
//	Feynman_Diagram_Ptr s_channel_3 = create_diagram("s_channel 3", Double_Wrench, VMP,
//	                                                 {Electron, Photon},
//	                                                 {Electron},
//	                                                 {Electron, Photon});
//
//	Feynman_Diagram_Ptr s_channel_4 = create_diagram("s_channel 4", Double_Wrench, VMP,
//	                                                 {Electron, Photon},
//	                                                 {Electron},
//	                                                 {Photon,Electron});


	Feynman_Diagram_Ptr t_channel = create_diagram("t_channel", X_Man, VMP,
	                                               {Electron, Photon},
	                                               {Electron},
	                                               {Photon, Electron});


	Feynman_Process e_scattering({t_channel});
//	Feynman_Process e_scattering_1({s_channel_1});
//	Feynman_Process e_scattering_2({s_channel_2});
//	Feynman_Process e_scattering_3({s_channel_3});
//	Feynman_Process e_scattering_4({s_channel_4});


	e_scattering.dsigma_dcos_table(std::cout, 500._MeV, {{0.2}});
//	e_scattering_1.dsigma_dcos_table(std::cout, 500._MeV, 0.1);
//	e_scattering_2.dsigma_dcos_table(std::cout, 500._MeV, 0.1);
//	e_scattering_3.dsigma_dcos_table(std::cout, 500._MeV, 0.1);
//	e_scattering_4.dsigma_dcos_table(std::cout, 500._MeV, 0.1);
    return EXIT_SUCCESS;
}

