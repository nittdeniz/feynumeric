#include <feynumeric/core.hpp>
#include <feynumeric/qed.hpp>
#include <feynumeric/particle_manager.hpp>

#include "effective_lagrangian_model.hpp"
#include "form_factors.hpp"

#include <iostream>
#include <iomanip>


int main(int argc, char** argv)
{
	using namespace Feynumeric;
	using namespace Feynumeric::Units;

	if( argc != 2 )
	{
		critical_error("Program expects filename as a single and only argument.");
	}
	Particle_Manager P((std::string(argv[1])));
	Particle_Ptr const& D1232pp  = P["D1232pp"];
	Particle_Ptr const& D1232p   = P["D1232p"];
	Particle_Ptr const& D1232n   = P["D1232n"];
	Particle_Ptr const& D1232m   = P["D1232m"];
	Particle_Ptr const& N1440p   = P.get("N1440p");
	Particle_Ptr const& N1440n   = P.get("N1440n");
	Particle_Ptr const& N1520p   = P.get("N1520p");
	Particle_Ptr const& N1520n   = P.get("N1520n");
	Particle_Ptr const& N1535p   = P["N1535p"];
	Particle_Ptr const& N1535n   = P["N1535n"];
	Particle_Ptr const& Proton   = P.get("proton");
	Particle_Ptr const& Neutron  = P.get("neutron");
	Particle_Ptr const& Pi_Plus  = P.get("pi+");
	Particle_Ptr const& Pi_Minus = P.get("pi-");
	Particle_Ptr const& Pi_Zero  = P.get("pi0");

	init_vertices(P);


	std::vector<Particle_Ptr> particles_Np = {N1440p, N1520p, N1535p};
	std::vector<Particle_Ptr> particles_Nn = {N1440n, N1520n, N1535n};

	for( std::size_t i = 0; i < particles_Np.size(); ++i )
	{
		auto& Np = particles_Np[i];
		auto& Nn = particles_Nn[i];
		Np->user_data("gRNpi", 1.);
		Np->user_data("form_factor", identity);
		Nn->user_data("gRNpi", 1.);
		Nn->user_data("form_factor", identity);
		auto channel_decay_Np1 = create_diagram(FORMAT("decay {} to proton pi0", Np->name()), Decay_1_to_2, VMP,
		                                        {Np},
		                                        {},
		                                        {Proton, Pi_Zero}
		);
		auto channel_decay_Np2 = create_diagram(FORMAT("decay {} to neutron pi+", Np->name()), Decay_1_to_2, VMP,
		                                        {Np},
		                                        {},
		                                        {Neutron, Pi_Plus}
		);


		auto channel_decay_Nn1 = create_diagram(FORMAT("decay {} to neutron pi0", Nn->name()), Decay_1_to_2, VMP,
		                                        {Nn},
		                                        {},
		                                        {Neutron, Pi_Zero}
		);
		auto channel_decay_Nn2 = create_diagram(FORMAT("decay {} to proton pi-", Nn->name()), Decay_1_to_2, VMP,
		                                        {Nn},
		                                        {},
		                                        {Proton, Pi_Minus}
		);

		Feynman_Process decay_Np1({channel_decay_Np1});
		Feynman_Process decay_Np2({channel_decay_Np2});

		Feynman_Process decay_Nn1({channel_decay_Nn1});
		Feynman_Process decay_Nn2({channel_decay_Nn2});

		auto w1 = decay_Np1.decay_width();
		auto w2 = decay_Np2.decay_width();

		auto w3 = decay_Nn1.decay_width();
		auto w4 = decay_Nn2.decay_width();


		double const literature_value1 = Np->width() * ( Np->user_data<double>("branching_N_pi_upper") +
		                                                    Np->user_data<double>("branching_N_pi_lower")) / 2.;
		double const literature_value2 = Nn->width() * ( Nn->user_data<double>("branching_N_pi_upper") +
		                                                    Nn->user_data<double>("branching_N_pi_lower")) / 2.;

		std::cout << FORMAT("g({} -> N + pi): ", Np->name()) << std::setw(10) << std::setprecision(10)<< std::sqrt(literature_value1 / ( w1 + w2 )) << "\n";
		std::cout << FORMAT("g({} -> N + pi): ", Nn->name()) << std::setw(10) << std::setprecision(10)<< std::sqrt(literature_value2 / ( w3 + w4 )) << "\n";
	}

	for( std::size_t i = 0; i < particles_Np.size(); ++i )
	{
		auto& Np = particles_Np[i];
		auto& Nn = particles_Nn[i];
		Np->user_data("gRNphoton", 1.);
		Nn->user_data("gRNphoton", 1.);

		auto channel_decay_Np1 = create_diagram(FORMAT("decay {} to proton photon", Np->name()), Decay_1_to_2, VMP,
		                                        {Np},
		                                        {},
		                                        {Proton, QED::Photon}
		);

		auto channel_decay_Nn1 = create_diagram(FORMAT("decay {} to neutron photon", Nn->name()), Decay_1_to_2, VMP,
		                                        {Nn},
		                                        {},
		                                        {Neutron, QED::Photon}
		);

		Feynman_Process decay_Np1({channel_decay_Np1});
		Feynman_Process decay_Nn1({channel_decay_Nn1});

		auto w1 = decay_Np1.decay_width();
		auto w2 = decay_Nn1.decay_width();


		double const literature_value1 = Np->width() * ( Np->user_data<double>("branching_proton_photon_upper") +
		                                                     Np->user_data<double>("branching_proton_photon_lower")) / 2.;
		double const literature_value2 = Nn->width() * ( Nn->user_data<double>("branching_neutron_photon_upper") +
		                                                     Nn->user_data<double>("branching_neutron_photon_lower")) / 2.;
		std::cout << FORMAT("g({} -> proton + photon): ", Np->name()) << std::setw(10) << std::setprecision(10)  << std::sqrt(literature_value1 / w1) << "\n";
		std::cout << FORMAT("g({} -> neutron + photon): ", Nn->name()) << std::setw(10) << std::setprecision(10) << std::sqrt(literature_value2 / w2) << "\n";
	}

	// coupling constant D1232 to N pi
	{
		D1232pp->user_data("gRNpi", 1.);
		D1232pp->user_data("form_factor", identity);
		D1232p->user_data("gRNpi", 1.);
		D1232p->user_data("form_factor", identity);
		D1232n->user_data("gRNpi", 1.);
		D1232n->user_data("form_factor", identity);
		D1232m->user_data("gRNpi", 1.);
		D1232m->user_data("form_factor", identity);
		auto channel_decay_D1232pp1 = create_diagram("decay_D1232++ to proton pi+", Decay_1_to_2, VMP,
		                                             {D1232pp},
		                                             {},
		                                             {Proton, Pi_Plus}
		);

		auto channel_decay_D1232p1 = create_diagram("decay_D1232+ to neutron pi+", Decay_1_to_2, VMP,
		                                            {D1232p},
		                                            {},
		                                            {Neutron, Pi_Plus}
		);
		auto channel_decay_D1232p2 = create_diagram("decay_ND1232+ to proton pi0", Decay_1_to_2, VMP,
		                                            {D1232p},
		                                            {},
		                                            {Proton, Pi_Zero}
		);

		auto channel_decay_D1232n1 = create_diagram("decay_D1232_0 to neutron pi0", Decay_1_to_2, VMP,
		                                            {D1232n},
		                                            {},
		                                            {Neutron, Pi_Zero}
		);
		auto channel_decay_D1232n2 = create_diagram("decay_ND1232_0 to proton pi-", Decay_1_to_2, VMP,
		                                            {D1232n},
		                                            {},
		                                            {Proton, Pi_Minus}
		);

		auto channel_decay_D1232m1 = create_diagram("decay_ND1232- to neutron pi-", Decay_1_to_2, VMP,
		                                            {D1232m},
		                                            {},
		                                            {Neutron, Pi_Minus}
		);

		Feynman_Process decay_D1232pp1({channel_decay_D1232pp1});

		Feynman_Process decay_D1232p1({channel_decay_D1232p1});
		Feynman_Process decay_D1232p2({channel_decay_D1232p2});

		Feynman_Process decay_D1232n1({channel_decay_D1232n1});
		Feynman_Process decay_D1232n2({channel_decay_D1232n2});

		Feynman_Process decay_D1232m1({channel_decay_D1232m1});

		auto w1 = decay_D1232pp1.decay_width();

		auto w2 = decay_D1232p1.decay_width();
		auto w3 = decay_D1232p2.decay_width();

		auto w4 = decay_D1232n1.decay_width();
		auto w5 = decay_D1232n2.decay_width();

		auto w6 = decay_D1232m1.decay_width();



		double const literature_value1 = D1232pp->width() * D1232pp->user_data<double>("branching_N_pi");
		double const literature_value2 = D1232p->width() * D1232p->user_data<double>("branching_N_pi");
		double const literature_value3 = D1232n->width() * D1232n->user_data<double>("branching_N_pi");
		double const literature_value4 = D1232m->width() * D1232m->user_data<double>("branching_N_pi");

		std::cout << "g(D1232pp -> N + pi): " << std::setw(10) << std::setprecision(10) << std::sqrt(literature_value1 / ( w1 )) << "\n";
		std::cout << "g(D1232p -> N + pi): "  << std::setw(10) << std::setprecision(10) << std::sqrt(literature_value2 / ( w2 + w3 )) << "\n";
		std::cout << "g(D1232n -> N + pi): "  << std::setw(10) << std::setprecision(10) << std::sqrt(literature_value3 / ( w4 + w5 )) << "\n";
		std::cout << "g(D1232m -> N + pi): "  << std::setw(10) << std::setprecision(10) << std::sqrt(literature_value4 / ( w6 )) << "\n";
	}

	// coupling constant D1232 to N gamma
	{
		D1232p->user_data("gRNphoton", 1.);
		D1232p->user_data("form_factor", identity);
		D1232n->user_data("gRNphoton", 1.);
		D1232n->user_data("form_factor", identity);

		auto channel_decay_D1232p1 = create_diagram("decay_D1232+ to proton photon", Decay_1_to_2, VMP,
		                                            {D1232p},
		                                            {},
		                                            {Proton, QED::Photon}
		);

		auto channel_decay_D1232n1 = create_diagram("decay_D1232_0 to neutron photon", Decay_1_to_2, VMP,
		                                            {D1232n},
		                                            {},
		                                            {Neutron, QED::Photon}
		);



		Feynman_Process decay_D1232p1({channel_decay_D1232p1});
		Feynman_Process decay_D1232n1({channel_decay_D1232n1});


		auto w1 = decay_D1232p1.decay_width();
		auto w2 = decay_D1232n1.decay_width();




		double const literature_value1 = D1232p->width() * ( D1232p->user_data<double>("branching_N_photon_upper") +
		                                                      D1232p->user_data<double>("branching_N_photon_lower")) / 2.;

		std::cout << "g(D1232p -> N + photon): "  << std::setw(10) << std::setprecision(10) << std::sqrt(literature_value1 / ( w1 )) << "\n";
		std::cout << "g(D1232n -> N + photon): "  << std::setw(10) << std::setprecision(10) << std::sqrt(literature_value1 / ( w2 )) << "\n";
	}
	return EXIT_SUCCESS;
}