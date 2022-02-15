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


    Feynman_Diagram_Ptr e_muon_t_channel = std::make_shared<Feynman_Diagram>(X_Man, VMP,
						   std::vector<Particle_Ptr>({Electron, Muon_Minus}),
						   std::vector<Particle_Ptr>({Photon}),
						   std::vector<Particle_Ptr>({Electron, Muon_Minus}));

	Feynman_Diagram_Ptr e_muon_t_channel2 = std::make_shared<Feynman_Diagram>(X_Man, VMP,
	                                                                         std::vector<Particle_Ptr>({Electron, Muon_Minus}),
	                                                                         std::vector<Particle_Ptr>({Photon}),
	                                                                         std::vector<Particle_Ptr>({Electron, Muon_Minus}));

	Feynman_Diagram_Ptr pair_production_s_channel = std::make_shared<Feynman_Diagram>(Double_Wrench, VMP,
                          std::vector<Particle_Ptr>({Photon, Positron}),
                            std::vector<Particle_Ptr>({Photon}),
                          std::vector<Particle_Ptr>({Muon_Plus, Muon_Minus}));


    Feynman_Process e_muon_scattering({1 * e_muon_t_channel, e_muon_t_channel2});
    e_muon_scattering.dsigma_dcos_table(std::cout, 1.49_GeV, 20ULL);
    return EXIT_SUCCESS;
}

