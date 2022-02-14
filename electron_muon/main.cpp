#include <iostream>


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


    Feynman_Diagram e_muon_t_channel(X_Man, VMP,
						   {Electron, Muon_Minus},
						   {Photon},
						   {Electron, Muon_Minus});

    Feynman_Diagram pair_production_s_channel(Double_Wrench, VMP,
                                 {Photon, Positron},
                                 {Photon},
                                 {Muon_Plus, Muon_Minus});


    Feynman_Process e_muon_scattering({&e_muon_t_channel});
    e_muon_scattering.dsigma_dcos_table(std::cout, 1.49_GeV, 20ULL);
    return EXIT_SUCCESS;
}

