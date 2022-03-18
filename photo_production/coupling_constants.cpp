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

	/********************************************\
	|*	███    ██  ██ ██   ██ ██   ██  ██████   *|
	|*	████   ██ ███ ██   ██ ██   ██ ██  ████  *|
	|*	██ ██  ██  ██ ███████ ███████ ██ ██ ██  *|
	|*	██  ██ ██  ██      ██      ██ ████  ██  *|
	|*	██   ████  ██      ██      ██  ██████   *|
    \********************************************/
	// coupling constant N1440 to N pi
	{
		N1440p->user_data("gRNpi", 1.);
		N1440p->user_data("gRNgamma", 1.);
		N1440p->user_data("gRNrho", 1.);
		N1440n->user_data("gRNpi", 1.);
		N1440n->user_data("gRNgamma", 1.);
		N1440n->user_data("gRNrho", 1.);
		auto channel_decay_n1440p1 = create_diagram("decay_n1440p to proton pi0", Wrench, VMP,
		                                            {N1440p},
		                                            {},
		                                            {Proton, Pi_Zero}
		);
		auto channel_decay_n1440p2 = create_diagram("decay_n1440p to neutron pi+", Wrench, VMP,
		                                            {N1440p},
		                                            {},
		                                            {Neutron, Pi_Plus}
		);

		auto channel_decay_n1440n1 = create_diagram("decay_n1440n to neutron pi0", Wrench, VMP,
		                                            {N1440n},
		                                            {},
		                                            {Neutron, Pi_Zero}
		);
		auto channel_decay_n1440n2 = create_diagram("decay_n1440n to proton pi-", Wrench, VMP,
		                                            {N1440n},
		                                            {},
		                                            {Proton, Pi_Minus}
		);

		Feynman_Process decay_n1440p1({channel_decay_n1440p1});
		Feynman_Process decay_n1440p2({channel_decay_n1440p2});

		Feynman_Process decay_n1440n1({channel_decay_n1440n1});
		Feynman_Process decay_n1440n2({channel_decay_n1440n2});

		auto w1 = decay_n1440p1.decay_width();
		auto w2 = decay_n1440p2.decay_width();

		auto w3 = decay_n1440n1.decay_width();
		auto w4 = decay_n1440n2.decay_width();


		double const literature_value1 = N1440p->width() * ( N1440p->user_data<long double>("branching_N_pi_upper") +
		                                                    N1440p->user_data<long double>("branching_N_pi_lower")) / 2.;
		double const literature_value2 = N1440n->width() * ( N1440n->user_data<long double>("branching_N_pi_upper") +
		                                                    N1440n->user_data<long double>("branching_N_pi_lower")) / 2.;

		std::cout << "g(N1440p -> N + pi): " << std::sqrt(literature_value1 / ( w1 + w2 )) << "\n";
		std::cout << "g(N1440n -> N + pi): " << std::sqrt(literature_value2 / ( w3 + w4 )) << "\n";
	}
	// coupling constant N1440 to N gamma
	{
		N1440p->user_data("gRNpi", 1.);
		N1440p->user_data("gRproton_photon", 1.);
		N1440p->user_data("gRNrho", 1.);
		N1440n->user_data("gRNpi", 1.);
		N1440n->user_data("gRneutron_photon", 1.);
		N1440n->user_data("gRNrho", 1.);

		auto channel_decay_n1440p1 = create_diagram("decay_n1440p to proton photon", Wrench, VMP,
		                                            {N1440p},
		                                            {},
		                                            {Proton, QED::Photon}
		);

		auto channel_decay_n1440n1 = create_diagram("decay_n1440n to neutron photon", Wrench, VMP,
		                                            {N1440n},
		                                            {},
		                                            {Neutron, QED::Photon}
		);

		//	channel_decay_n1440p->print_amplitude();
		//
		Feynman_Process decay_n1440p1({channel_decay_n1440p1});
		Feynman_Process decay_n1440n1({channel_decay_n1440n1});

		auto w1 = decay_n1440p1.decay_width();
		auto w2 = decay_n1440n1.decay_width();


		double const literature_value1 = N1440p->width() * ( N1440p->user_data<long double>("branching_proton_photon_upper") +
		                                                     N1440p->user_data<long double>("branching_proton_photon_lower")) / 2.;
		double const literature_value2 = N1440n->width() * ( N1440n->user_data<long double>("branching_neutron_photon_upper") +
		                                                     N1440n->user_data<long double>("branching_neutron_photon_lower")) / 2.;
		std::cout << "g(N1440p -> proton + photon): " << std::sqrt(literature_value1 / w1) << "\n";
		std::cout << "g(N1440n -> neutron + photon): " << std::sqrt(literature_value2 / w2) << "\n";
	}

	/********************************************\
	|*	███    ██  ██ ███████ ██████   ██████   *|
	|*	████   ██ ███ ██           ██ ██  ████  *|
	|*	██ ██  ██  ██ ███████  █████  ██ ██ ██  *|
	|*	██  ██ ██  ██      ██ ██      ████  ██  *|
	|*	██   ████  ██ ███████ ███████  ██████   *|
	\********************************************/
	// coupling constant N1520 to N pi
	{
		N1520p->user_data("gRNpi", 1.);
		N1520p->user_data("gRNgamma", 1.);
		N1520p->user_data("gRNrho", 1.);
		N1520n->user_data("gRNpi", 1.);
		N1520n->user_data("gRNgamma", 1.);
		N1520n->user_data("gRNrho", 1.);
		auto channel_decay_N1520p1 = create_diagram("decay_N1520p to proton pi0", Wrench, VMP,
		                                            {N1520p},
		                                            {},
		                                            {Proton, Pi_Zero}
		);
		auto channel_decay_N1520p2 = create_diagram("decay_N1520p to neutron pi+", Wrench, VMP,
		                                            {N1520p},
		                                            {},
		                                            {Neutron, Pi_Plus}
		);

		auto channel_decay_N1520n1 = create_diagram("decay_N1520n to neutron pi0", Wrench, VMP,
		                                            {N1520n},
		                                            {},
		                                            {Neutron, Pi_Zero}
		);
		auto channel_decay_N1520n2 = create_diagram("decay_N1520n to proton pi-", Wrench, VMP,
		                                            {N1520n},
		                                            {},
		                                            {Proton, Pi_Minus}
		);

		Feynman_Process decay_N1520p1({channel_decay_N1520p1});
		Feynman_Process decay_N1520p2({channel_decay_N1520p2});

		Feynman_Process decay_N1520n1({channel_decay_N1520n1});
		Feynman_Process decay_N1520n2({channel_decay_N1520n2});

		auto w1 = decay_N1520p1.decay_width();
		auto w2 = decay_N1520p2.decay_width();

		auto w3 = decay_N1520n1.decay_width();
		auto w4 = decay_N1520n2.decay_width();


		double const literature_value1 = N1520p->width() * ( N1520p->user_data<long double>("branching_N_pi_upper") +
		                                                     N1520p->user_data<long double>("branching_N_pi_lower")) / 2.;
		double const literature_value2 = N1520n->width() * ( N1520n->user_data<long double>("branching_N_pi_upper") +
		                                                     N1520n->user_data<long double>("branching_N_pi_lower")) / 2.;


		std::cout << "g(N1520p -> N + pi): " << std::sqrt(literature_value1 / ( w1 + w2 )) << "\n";
		std::cout << "g(N1520n -> N + pi): " << std::sqrt(literature_value2 / ( w3 + w4 )) << "\n";
	}
	/*
	// coupling constant N1520 to N gamma
	{
		N1440p->user_data("gRNpi", 1.);
		N1440p->user_data("gRproton_photon", 1.);
		N1440p->user_data("gRNrho", 1.);
		N1440n->user_data("gRNpi", 1.);
		N1440n->user_data("gRneutron_photon", 1.);
		N1440n->user_data("gRNrho", 1.);

		auto channel_decay_n1440p1 = create_diagram("decay_n1440p to proton photon", Wrench, VMP,
		                                            {N1440p},
		                                            {},
		                                            {Proton, QED::Photon}
		);

		auto channel_decay_n1440n1 = create_diagram("decay_n1440n to neutron photon", Wrench, VMP,
		                                            {N1440n},
		                                            {},
		                                            {Neutron, QED::Photon}
		);

		//	channel_decay_n1440p->print_amplitude();
		//
		Feynman_Process decay_n1440p1({channel_decay_n1440p1});
		Feynman_Process decay_n1440n1({channel_decay_n1440n1});

		auto w1 = decay_n1440p1.decay_width();
		auto w2 = decay_n1440n1.decay_width();


		double const literature_value1 = N1440p->width() * ( N1440p->user_data<long double>("branching_proton_photon_upper") +
		                                                     N1440p->user_data<long double>("branching_proton_photon_lower")) / 2.;
		double const literature_value2 = N1440n->width() * ( N1440n->user_data<long double>("branching_neutron_photon_upper") +
		                                                     N1440n->user_data<long double>("branching_neutron_photon_lower")) / 2.;
		std::cout << "g(N1520p -> proton + photon): " << std::sqrt(literature_value1 / w1) << "\n";
		std::cout << "g(N1520n -> neutron + photon): " << std::sqrt(literature_value2 / w2) << "\n";
	}
	 */
	return EXIT_SUCCESS;
}