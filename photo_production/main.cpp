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

	N1440p->user_data("gRNPi", 0.35*0.65 / decay_n1440p.decay_width());

	return EXIT_SUCCESS;
}