#include <feynumeric/core.hpp>
#include <feynumeric/qed.hpp>

#include "effective_lagrangian_model.hpp"

#include <iomanip>
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
		auto channel_s1 = create_diagram("N1440 s_channel1", Double_Wrench, VMP,
		                                {Proton, QED::Photon},
		                                {N1440p},
		                                {Pi_Zero, Proton}
		);
		auto channel_u1 = create_diagram("N1440 u_channel1", X_Man, VMP,
		                                 {Proton, QED::Photon},
		                                 {N1440p},
		                                 {Pi_Zero, Proton}
		);

		auto channel_s2 = create_diagram("N1440 s_channel2", Double_Wrench, VMP,
		                                {Neutron, QED::Photon},
		                                {N1440n},
		                                {Pi_Minus, Proton}
		);

		auto channel_u2 = create_diagram("N1440 u_channel2", X_Man, VMP,
		                                 {Neutron, QED::Photon},
		                                 {N1440p},
		                                 {Pi_Minus, Proton}
		);

		auto channel_s3 = create_diagram("N1440 s_channel3", Double_Wrench, VMP,
		                                {Proton, QED::Photon},
		                                {N1440p},
		                                {Pi_Plus, Neutron}
		);

		auto channel_u3 = create_diagram("N1440 u_channel3", X_Man, VMP,
		                                 {Proton, QED::Photon},
		                                 {N1440n},
		                                 {Pi_Plus, Neutron}
		);


		Feynman_Process scattering_n1440p1({channel_s1, channel_u1});
		Feynman_Process scattering_n1440p2({channel_s2, channel_u2});
		Feynman_Process scattering_n1440p3({channel_s3, channel_u3});

//		std::cout << "1.2GeV\n";
//
//		scattering_n1440p1.print_dsigma_dcos_table(std::cout, 1.2_GeV, 0.1);
//		scattering_n1440p2.print_dsigma_dcos_table(std::cout, 1.2_GeV, 0.1);
//		scattering_n1440p3.print_dsigma_dcos_table(std::cout, 1.2_GeV, 0.1);
//
//		std::cout << "1.6GeV\n";
//
//		scattering_n1440p1.print_dsigma_dcos_table(std::cout, 1.6_GeV, 0.1);
//		scattering_n1440p2.print_dsigma_dcos_table(std::cout, 1.6_GeV, 0.1);
//		scattering_n1440p3.print_dsigma_dcos_table(std::cout, 1.6_GeV, 0.1);

		std::cout << "2GeV\n";
		scattering_n1440p1.print_dsigma_dcos_table(std::cout, 2._GeV, 0.1);
		scattering_n1440p2.print_dsigma_dcos_table(std::cout, 2._GeV, 0.1);
		scattering_n1440p3.print_dsigma_dcos_table(std::cout, 2._GeV, 0.1);
		//scattering_n1440p.print_sigma_table(std::cout, {{1.5_GeV}});
	}
	return EXIT_SUCCESS;
}