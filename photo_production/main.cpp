#include <feynumeric/core.hpp>
#include <feynumeric/qed.hpp>

#include "effective_lagrangian_model.hpp"

#include <iostream>

int main()
{
	using namespace Feynumeric;
	using namespace Feynumeric::Units;

	init_particles();
	init_vertices();

	// coupling constant N1440p to Proton pi0
	auto channel_decay_n1440p = create_diagram("decay_n1440p to proton pi0", Wrench, VMP,
	                                   {N1440p},
	                                   {},
	                                   {Proton, Pi_Zero}
									);
	Feynman_Process decay_n1440p({channel_decay_n1440p});

	std::cout << "g: " << std::sqrt(0.35*0.65 / decay_n1440p.decay_width()) << "\n";

	N1440p->user_data("gRNpi", 0.38);
	N1440p->user_data("gRNgamma", 0.0528755);

	// Photo Production
	{
		auto channel_s = create_diagram("N1440 s_channel", Double_Wrench, VMP,
		                                {Proton, QED::Photon},
		                                {N1440p},
		                                {Pi_Zero, Proton}
		);

		Feynman_Process scattering_n1440p({channel_s});
		scattering_n1440p.dsigma_dcos_table(std::cout, 1.49_GeV, 0.1);
		scattering_n1440p.sigma_table(std::cout, {{1.4_GeV}});
	}

	return EXIT_SUCCESS;
}