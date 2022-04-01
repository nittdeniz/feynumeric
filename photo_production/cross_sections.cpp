#include <feynumeric/core.hpp>
#include <feynumeric/qed.hpp>

#include "effective_lagrangian_model.hpp"

#include <iostream>

int main(int argc, char** argv)
{
	using namespace Feynumeric;
	using namespace Feynumeric::Units;

	if( argc != 2 )
	{
		critical_error("Program expects filename as a single and only argument.");
	}
	Particle_Manager P((std::string(argv[1])));
	Particle_Ptr const& N1440p   = P.get("N1440p");
	Particle_Ptr const& N1440n   = P.get("N1440n");
	Particle_Ptr const& N1520p   = P.get("N1520p");
	Particle_Ptr const& N1520n   = P.get("N1520n");
	Particle_Ptr const& Proton   = P.get("proton");
	Particle_Ptr const& Neutron  = P.get("neutron");
	Particle_Ptr const& Pi_Plus  = P.get("pi+");
	Particle_Ptr const& Pi_Minus = P.get("pi-");
	Particle_Ptr const& Pi_Zero  = P.get("pi0");

	init_vertices(P);

	//init_particles();

//	N1440p->user_data("gRNpi", 0.38);
//	N1440p->user_data("gRNgamma", 0.0528755);

	// Photo Production
	{
		auto channel_s = create_diagram("N1440 s_channel", Double_Wrench, VMP,
		                                {Proton, QED::Photon},
		                                {N1440p},
		                                {Pi_Zero, Proton}
		);


		Feynman_Process scattering_n1440p({channel_s});
		scattering_n1440p.dsigma_dcos_table(std::cout, 1.49_GeV, 0.1);
//		scattering_n1440p.sigma_table(std::cout, {{1.4_GeV}});
	}
	return EXIT_SUCCESS;
}