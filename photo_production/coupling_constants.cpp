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
	auto const& Proton = P["proton"];
	auto const& Neutron = P["neutron"];
	auto const& Photon = P["photon"];
	auto const& Pi_Zero = P["pi0"];
	auto const& Pi_Minus = P["pi-"];
	auto const& Pi_Plus = P["pi+"];

	init_vertices(P);

	std::vector<std::string> particle_strings = {"N1440", "N1520", "N1535", "N1650", "N1675", "N1680", "N1700", "N1710", "N1720", "N1875", "N1880", "N1895", "N1900", "N2060", "N2100", "N2120","D1232", "D1600", "D1620", "D1700", "D1900", "D1910", "D1920"};

	std::vector<Particle_Ptr> particles_Np;
	std::vector<Particle_Ptr> particles_Nn;
	std::vector<Particle_Ptr> particles_Dpp;
	std::vector<Particle_Ptr> particles_Dm;

	for( auto const& str : particle_strings ){
		if( str[0] == 'D' ){
			particles_Dpp.push_back(P.get(str + "pp"));
			particles_Dm.push_back(P.get(str + "m"));
		}
		particles_Np.push_back(P.get(str + "p"));
		particles_Nn.push_back(P.get(str + "n"));
	}


//	std::vector<Particle_Ptr> particles_Np = {P["N1440p"], P["N1520p"], P["N1535p"],P["D1232p"],  P["D1600p"], P["D1920p"]};
//	std::vector<Particle_Ptr> particles_Nn = {P["N1440n"], P["N1520n"], P["N1535n"],P["D1232n"],  P["D1600n"], P["D1920n"]};

//	std::vector<Particle_Ptr> particles_Dpp = {P["D1232pp"], P["D1600pp"], P["D1920pp"]};
//	std::vector<Particle_Ptr> particles_Dp  = {P["D1232p"],  P["D1600p"], P["D1920p"]};
//	std::vector<Particle_Ptr> particles_Dn  = {P["D1232n"],  P["D1600n"], P["D1920n"]};
//	std::vector<Particle_Ptr> particles_Dm  = {P["D1232m"],  P["D1600m"], P["D1920m"]};



	for( std::size_t i = 0; i < particles_Np.size(); ++i )
	{
		auto& Np = particles_Np[i];
		auto& Nn = particles_Nn[i];
		if( Np->spin().j() > 2 ){
			continue;
		}
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
		if( Np->spin().j() > 2 ){
			continue;
		}
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

	for( std::size_t i = 0; i < particles_Dpp.size(); ++i )
	{
		auto& Dpp = particles_Dpp[i];
//		auto& Dp = particles_Dp[i];
//		auto& Dn = particles_Dn[i];
		auto& Dm = particles_Dm[i];
		if( Dpp->spin().j() > 2 ){
			continue;
		}
		Dpp->user_data("gRNpi", 1.);
		Dpp->user_data("form_factor", identity);
//		Dp->user_data("gRNpi", 1.);
//		Dp->user_data("form_factor", identity);
//		Dn->user_data("gRNpi", 1.);
//		Dn->user_data("form_factor", identity);
		Dm->user_data("gRNpi", 1.);
		Dm->user_data("form_factor", identity);

		auto channel_decay_Dpp1 = create_diagram(FORMAT("decay {} to proton pi+", Dpp->name()), Decay_1_to_2, VMP,
		                                             {Dpp},
		                                             {},
		                                             {Proton, Pi_Plus}
		);

//		auto channel_decay_Dp1 = create_diagram(FORMAT("decay {} to neutron pi+", Dp->name()), Decay_1_to_2, VMP,
//		                                            {Dp},
//		                                            {},
//		                                            {Neutron, Pi_Plus}
//		);
//		auto channel_decay_Dp2 = create_diagram(FORMAT("decay {} to proton pi0", Dp->name()), Decay_1_to_2, VMP,
//		                                            {Dp},
//		                                            {},
//		                                            {Proton, Pi_Zero}
//		);
//
//		auto channel_decay_Dn1 = create_diagram(FORMAT("decay {} to neutron pi0", Dn->name()), Decay_1_to_2, VMP,
//		                                            {Dn},
//		                                            {},
//		                                            {Neutron, Pi_Zero}
//		);
//		auto channel_decay_Dn2 = create_diagram(FORMAT("decay {} to proton pi-", Dn->name()), Decay_1_to_2, VMP,
//		                                            {Dn},
//		                                            {},
//		                                            {Proton, Pi_Minus}
//		);

		auto channel_decay_Dm1 = create_diagram(FORMAT("decay {} to neutron pi-", Dm->name()), Decay_1_to_2, VMP,
		                                            {Dm},
		                                            {},
		                                            {Neutron, Pi_Minus}
		);

		Feynman_Process decay_Dpp1({channel_decay_Dpp1});
		Feynman_Process decay_Dm1({channel_decay_Dm1});

		auto w1 = decay_Dpp1.decay_width();
		auto w6 = decay_Dm1.decay_width();
		double const literature_value1 = Dpp->width() *  ( Dpp->user_data<double>("branching_N_pi_upper") +
		                                                   Dpp->user_data<double>("branching_N_pi_lower")) / 2.;
		double const literature_value4 = Dm->width() *  ( Dm->user_data<double>("branching_N_pi_upper") +
		                                                  Dm->user_data<double>("branching_N_pi_lower")) / 2.;

		std::cout << FORMAT("g({} -> N + pi): ", Dpp->name()) << std::setw(10) << std::setprecision(10) << std::sqrt(literature_value1 / ( w1 )) << "\n";
//		std::cout << FORMAT("g({} -> N + pi): ", Dp->name())  << std::setw(10) << std::setprecision(10) << std::sqrt(literature_value2 / ( w2 + w3 )) << "\n";
//		std::cout << FORMAT("g({} -> N + pi): ", Dn->name())  << std::setw(10) << std::setprecision(10) << std::sqrt(literature_value3 / ( w4 + w5 )) << "\n";
		std::cout << FORMAT("g({} -> N + pi): ", Dm->name())  << std::setw(10) << std::setprecision(10) << std::sqrt(literature_value4 / ( w6 )) << "\n";
	}

//	for( std::size_t i = 0; i < particles_Dpp.size(); ++i )
//	{
//		auto& Dpp = particles_Dpp[i];
////		auto& Dp = particles_Dp[i];
////		auto& Dn = particles_Dn[i];
//		auto& Dm = particles_Dm[i];
//		Dpp->user_data("gRNphoton", 1.);
//		Dpp->user_data("form_factor", identity);
////		Dp->user_data("gRNphoton", 1.);
////		Dp->user_data("form_factor", identity);
////		Dn->user_data("gRNphoton", 1.);
////		Dn->user_data("form_factor", identity);
//		Dm->user_data("gRNphoton", 1.);
//		Dm->user_data("form_factor", identity);
//
//		auto channel_decay_Dp1 = create_diagram(FORMAT("decay {} to proton photon", Dp->name()), Decay_1_to_2, VMP,
//		                                            {Dp},
//		                                            {},
//		                                            {Proton, QED::Photon}
//		);
//
//		auto channel_decay_Dn1 = create_diagram(FORMAT("decay {} to neutron photon", Dn->name()), Decay_1_to_2, VMP,
//		                                            {Dn},
//		                                            {},
//		                                            {Neutron, QED::Photon}
//		);
//
//
//
//		Feynman_Process decay_Dp1({channel_decay_Dp1});
//		Feynman_Process decay_Dn1({channel_decay_Dn1});
//
//
//		auto w1 = decay_Dp1.decay_width();
//		auto w2 = decay_Dn1.decay_width();
//
//
//
//
//		double const literature_value1 = Dp->width() * ( Dp->user_data<double>("branching_N_photon_upper") +
//		                                                      Dp->user_data<double>("branching_N_photon_lower")) / 2.;
//
//		std::cout << FORMAT("g({} -> N + photon): ", Dp->name())  << std::setw(10) << std::setprecision(10) << std::sqrt(literature_value1 / ( w1 )) << "\n";
//		std::cout << FORMAT("g({} -> N + photon): ", Dn->name())  << std::setw(10) << std::setprecision(10) << std::sqrt(literature_value1 / ( w2 )) << "\n";
//	}
	return EXIT_SUCCESS;
}