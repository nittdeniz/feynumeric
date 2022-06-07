#include <feynumeric/core.hpp>
#include <feynumeric/qed.hpp>
#include <feynumeric/particle_manager.hpp>

#include "effective_lagrangian_model.hpp"
#include "form_factors.hpp"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

#define Mesons 1
#define Npi 1
#define Nphoton 1
#define N1440 1
#define N1520 1
#define N1535 1
#define D1600 1
#define D1620 1
#define D1700 1


int main(int argc, char** argv)
{
	using namespace Feynumeric;
	using namespace Feynumeric::Units;

	if( argc != 2 )
	{
		critical_error("Program expects filename as a single and only argument.");
	}

	std::stringstream buffer;

	Particle_Manager P((std::string(argv[1])));
	auto const& Proton = P.get("proton");
	auto const& Neutron = P.get("neutron");
	auto const& Pi_Zero = P.get("pi0");
	auto const& Pi_Minus = P.get("pi-");
	auto const& Pi_Plus = P.get("pi+");

	for( auto& [key, particle] : P ){
		particle->user_data("form_factor", identity);
	}

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
#if Npi
    for( std::size_t i = 0; i < particles_Np.size(); ++i )
    {
        auto& Np = particles_Np[i];
        if( Np->spin().j() > 2 ){
            continue;
        }
        auto coupl_str = coupling_string(Np->name().substr(0, 5), "N", "Pion");
        couplings.set(coupl_str, 1.);

        if( !Np->exists_user_data("branching_N_pi") ){
            Np->user_data("branching_N_pi", 1.);
        }
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

        Feynman_Process decay_Np1({channel_decay_Np1});
        Feynman_Process decay_Np2({channel_decay_Np2});

        auto w1 = decay_Np1.decay_width();
        auto w2 = decay_Np2.decay_width();

        auto const literature_value = Np->width() * Np->user_data<double>("branching_N_pi");
        auto const g = std::sqrt(literature_value / ( w1 + w2 ));
        couplings.set(coupl_str, g);
        buffer <<  FORMAT("{} {}\n", coupl_str, g);
    }
#endif
#if Nphoton
    for( std::size_t i = 0; i < particles_Np.size(); ++i )
    {
        auto& Np = particles_Np[i];
        auto& Nn = particles_Nn[i];
        if( Np->spin().j() > 2 ){
            continue;
        }
        auto coupl_str1 = coupling_string(Np->name().substr(0, 5), Proton->name(), QED::Photon->name());
        auto coupl_str2 = coupling_string(Nn->name().substr(0, 5), Neutron->name(), QED::Photon->name());
        couplings.set(coupl_str1, 1.);
        couplings.set(coupl_str2, 1.);

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


        double const literature_value1 = Np->width() * Np->user_data<double>("branching_proton_photon");
        double const literature_value2 = Nn->width() * Nn->user_data<double>("branching_neutron_photon");
        auto g1 = std::sqrt(literature_value1 / w1);
        auto g2 = std::sqrt(literature_value2 / w2);
        couplings.set(coupl_str1, g1);
        couplings.set(coupl_str2, g2);
        buffer <<  FORMAT("{} {}\n", coupl_str1, g1);
        buffer <<  FORMAT("{} {}\n", coupl_str2, g2);
    }
#endif

#if Npi
    for( std::size_t i = 0; i < particles_Dpp.size(); ++i )
    {
        auto& Dpp = particles_Dpp[i];
        if( Dpp->spin().j() > 2 ){
            continue;
        }
        auto coupl_str = coupling_string(Dpp->name().substr(0, 5), "N", "Pion");
        couplings.set(coupl_str, 1.);

        auto channel_decay_Dpp1 = create_diagram(FORMAT("decay {} to proton pi+", Dpp->name()), Decay_1_to_2, VMP,
                                                 {Dpp},
                                                 {},
                                                 {Proton, Pi_Plus}
        );
        Feynman_Process decay_Dpp1({channel_decay_Dpp1});
        auto w1 = decay_Dpp1.decay_width();
        auto const literature_value = Dpp->width() *  Dpp->user_data<double>("branching_N_pi");
        auto const g = std::sqrt(literature_value / ( w1 ));
        couplings.set(coupl_str, g);
        buffer <<  FORMAT("{} {}\n", coupl_str, g);
    }
#endif
	/* Other */
    #if Mesons
    {   /// f0(500) to pi pi
        auto particle = P.get("f0_500");
        auto coupl_str = coupling_string("f0_500", "Pion", "Pion");
        couplings.set(coupl_str, 1.);
        auto channel_decay_f0_1 = create_diagram(FORMAT("decay {} to pi+ pi-", particle->name()), Decay_1_to_2, VMP,
                                                 {particle},
                                                 {},
                                                 {Pi_Plus, Pi_Minus}
        );
        auto channel_decay_f0_2 = create_diagram(FORMAT("decay {} to pi0 pi0", particle->name()), Decay_1_to_2, VMP,
                                                 {particle},
                                                 {},
                                                 {Pi_Zero, Pi_Zero}
        );

        Feynman_Process decay1({channel_decay_f0_1});
        Feynman_Process decay2({channel_decay_f0_2});
        auto g = std::sqrt(particle->width() / (decay1.decay_width() + decay2.decay_width()));
        couplings.set(coupl_str, g);
        buffer <<  FORMAT("{} {}\n", coupl_str, g);
    }
    {   /// rho to pi pi
        auto particle = P.get("rho0");
        auto coupl_str = coupling_string("Rho", "Pion", "Pion");
        couplings.set(coupl_str, 1.);
        auto channel_decay_rho_1 = create_diagram(FORMAT("decay {} to pi+ pi-", particle->name()), Decay_1_to_2, VMP,
                                                  {particle},
                                                  {},
                                                  {Pi_Plus, Pi_Minus}
        );
        Feynman_Process decay1({channel_decay_rho_1});
        auto g = std::sqrt(particle->width() / (decay1.decay_width() ));
        couplings.set(coupl_str, g);
        buffer << FORMAT("{} {}\n", coupl_str, g);
    }
    #endif
    #if N1440
    {   /// N1440 to D1232
        auto coupl_str = coupling_string("N1440", "D1232", "Pion");
        couplings.set(coupl_str, 1.);

        auto channel_decay_N1440p_D1232_pi_1 = create_diagram(FORMAT("decay N1440 to Delta pi 1"), Decay_1_to_M2_1, VMP,
                                                         {P.get("N1440p")},
                                                         {P.get("D1232pp")},
                                                         {Proton, Pi_Plus, Pi_Minus}
        );
        auto channel_decay_N1440p_D1232_pi_2 = create_diagram(FORMAT("decay N1440 to Delta pi 1"), Decay_1_to_M2_1_cross, VMP,
                                                         {P.get("N1440p")},
                                                         {P.get("D1232n")},
                                                         {Proton, Pi_Plus, Pi_Minus}
        );

        auto channel_decay_N1440p_D1232_pi_3 = create_diagram(FORMAT("decay N1440 to Delta pi 1"), Decay_1_to_M2_1, VMP,
                                                              {P.get("N1440p")},
                                                              {P.get("D1232p")},
                                                              {Proton, Pi_Zero, Pi_Zero}
        );
        auto channel_decay_N1440p_D1232_pi_4 = create_diagram(FORMAT("decay N1440 to Delta pi 1"), Decay_1_to_M2_1_cross, VMP,
                                                              {P.get("N1440p")},
                                                              {P.get("D1232p")},
                                                              {Proton, Pi_Zero, Pi_Zero}
        );

        auto channel_decay_N1440p_D1232_pi_5 = create_diagram(FORMAT("decay N1440 to Delta pi 1"), Decay_1_to_M2_1, VMP,
                                                              {P.get("N1440p")},
                                                              {P.get("D1232p")},
                                                              {Neutron, Pi_Plus, Pi_Zero}
        );
        auto channel_decay_N1440p_D1232_pi_6 = create_diagram(FORMAT("decay N1440 to Delta pi 1"), Decay_1_to_M2_1_cross, VMP,
                                                              {P.get("N1440p")},
                                                              {P.get("D1232n")},
                                                              {Neutron, Pi_Plus, Pi_Zero}
        );

        auto channel_decay_N1440n_D1232_pi_5 = create_diagram(FORMAT("decay N1440 to Delta pi 1"), Decay_1_to_M2_1, VMP,
                                                              {P.get("N1440n")},
                                                              {P.get("D1232n")},
                                                              {Neutron, Pi_Zero, Pi_Zero}
        );
        auto channel_decay_N1440n_D1232_pi_6 = create_diagram(FORMAT("decay N1440 to Delta pi 1"), Decay_1_to_M2_1_cross, VMP,
                                                              {P.get("N1440n")},
                                                              {P.get("D1232n")},
                                                              {Neutron, Pi_Zero, Pi_Zero}
        );

        Feynman_Process decay_Np1({channel_decay_N1440p_D1232_pi_1, channel_decay_N1440p_D1232_pi_2});
        Feynman_Process decay_Np2({channel_decay_N1440p_D1232_pi_3, channel_decay_N1440p_D1232_pi_4});
        Feynman_Process decay_Np3({channel_decay_N1440p_D1232_pi_5, channel_decay_N1440p_D1232_pi_6});

        auto w1 = decay_Np1.decay_width();
        auto w2 = decay_Np2.decay_width();
        auto w3 = decay_Np3.decay_width();

        double const literature_value = P.get("N1440")->width() * P.get("N1440")->user_data<double>("branching_N_pipi_D1232");

        auto const g = std::sqrt(literature_value / ( w1 + w2 + w3 ));
        couplings.set(coupl_str, g);
        buffer << FORMAT("{} {}\n", coupl_str, g);
    }
	{   /// N1440 to f0(500)
        auto coupl_str = coupling_string("N1440", "N", "f0_500");
        couplings.set(coupl_str, 1.);
		auto channel_decay_N1440p_f500_1 = create_diagram(FORMAT("decay N1440 to f_0(500) pi+pi-"), Decay_1_to_M2_1, VMP,
		                                                      {P.get("N1440p")},
		                                                      {P.get("f0_500")},
		                                                      {Pi_Plus, Pi_Minus, Proton}
		);
		auto channel_decay_N1440p_f500_2 = create_diagram(FORMAT("decay N1440 to f_0(500) pi0pi0"), Decay_1_to_M2_1, VMP,
		                                                  {P.get("N1440p")},
		                                                  {P.get("f0_500")},
		                                                  {Pi_Zero, Pi_Zero, Proton}
		);

		auto const literature_value = P.get("N1440")->width() * P.get("N1440")->user_data<double>("branching_N_pipi_f500");

		Feynman_Process decay_Np1({channel_decay_N1440p_f500_1});
		Feynman_Process decay_Np2({channel_decay_N1440p_f500_2});

		auto const w1 = decay_Np1.decay_width();
		auto const w2 = decay_Np2.decay_width();

		auto const g = std::sqrt(literature_value / ( w1 + w2 ));
        couplings.set(coupl_str, g);
		buffer << FORMAT("{} {}\n", coupl_str, g);
	}
    #endif
    #if N1520
	{   /// N1520 to rho
        auto coupl_str = coupling_string("N1520", "N", "Rho");
        couplings.set(coupl_str, 1.);
		auto channel_decay_N1520p_rho_1 = create_diagram(FORMAT("decay N1520p to rho0 pi+pi-"), Decay_1_to_M2_1, VMP,
		                                                  {P.get("N1520p")},
		                                                  {P.get("rho0")},
		                                                  {Pi_Plus, Pi_Minus, Proton}
		);
		auto channel_decay_N1520p_rho_2 = create_diagram(FORMAT("decay N1520p to rho0 pi0pi0"), Decay_1_to_M2_1, VMP,
		                                                 {P.get("N1520p")},
		                                                 {P.get("rho+")},
		                                                 {Pi_Plus, Pi_Zero, Neutron}
		);
		double const literature_value = P.get("N1520")->width() * P.get("N1520")->user_data<double>("branching_N_pipi_rho");

		Feynman_Process decay_Np1({channel_decay_N1520p_rho_1});
		Feynman_Process decay_Np2({channel_decay_N1520p_rho_2});
		auto w1 = decay_Np1.decay_width();
		auto w2 = decay_Np2.decay_width();

        auto const g = std::sqrt(literature_value / ( w1 + w2 ));
        couplings.set(coupl_str, g);
        buffer << FORMAT("{} {}\n", coupl_str, g);
	}
	{   /// N1520 to D1232
        auto coupl_str = coupling_string("N1520", "D1232", "Pion");
        couplings.set(coupl_str, 1.);

		auto channel_decay_N1520p_D1232_pi_1 = create_diagram(FORMAT("decay N1520 to Delta pi 1"), Decay_1_to_M2_1, VMP,
		                                                      {P.get("N1520p")},
		                                                      {P.get("D1232pp")},
		                                                      {Proton, Pi_Plus, Pi_Minus}
		);
		auto channel_decay_N1520p_D1232_pi_2 = create_diagram(FORMAT("decay N1520 to Delta pi 1"), Decay_1_to_M2_1_cross, VMP,
		                                                      {P.get("N1520p")},
		                                                      {P.get("D1232n")},
		                                                      {Proton, Pi_Plus, Pi_Minus}
		);

		auto channel_decay_N1520p_D1232_pi_3 = create_diagram(FORMAT("decay N1520 to Delta pi 1"), Decay_1_to_M2_1, VMP,
		                                                      {P.get("N1520p")},
		                                                      {P.get("D1232p")},
		                                                      {Proton, Pi_Zero, Pi_Zero}
		);
		auto channel_decay_N1520p_D1232_pi_4 = create_diagram(FORMAT("decay N1520 to Delta pi 1"), Decay_1_to_M2_1_cross, VMP,
		                                                      {P.get("N1520p")},
		                                                      {P.get("D1232p")},
		                                                      {Proton, Pi_Zero, Pi_Zero}
		);

		auto channel_decay_N1520p_D1232_pi_5 = create_diagram(FORMAT("decay N1520 to Delta pi 1"), Decay_1_to_M2_1, VMP,
		                                                      {P.get("N1520p")},
		                                                      {P.get("D1232p")},
		                                                      {Neutron, Pi_Plus, Pi_Zero}
		);
		auto channel_decay_N1520p_D1232_pi_6 = create_diagram(FORMAT("decay N1520 to Delta pi 1"), Decay_1_to_M2_1_cross, VMP,
		                                                      {P.get("N1520p")},
		                                                      {P.get("D1232n")},
		                                                      {Neutron, Pi_Plus, Pi_Zero}
		);

		Feynman_Process decay_Np1({channel_decay_N1520p_D1232_pi_1, channel_decay_N1520p_D1232_pi_2});
		Feynman_Process decay_Np2({channel_decay_N1520p_D1232_pi_3, channel_decay_N1520p_D1232_pi_4});
		Feynman_Process decay_Np3({channel_decay_N1520p_D1232_pi_5, channel_decay_N1520p_D1232_pi_6});

		auto w1 = decay_Np1.decay_width();
		auto w2 = decay_Np2.decay_width();
		auto w3 = decay_Np3.decay_width();

		double const literature_value = P.get("N1520")->width() * P.get("N1520")->user_data<double>("branching_N_pipi");

        auto const g = std::sqrt(literature_value / ( w1 + w2 + w3 ));
        couplings.set(coupl_str, g);
        buffer << FORMAT("{} {}\n", coupl_str, g);
	}
    #endif
    #if N1535
    {   /// N1535 to eta
        auto coupl_str = coupling_string("N1535", "N", "eta");
        couplings.set(coupl_str, 1.);
        auto channel_decay_N1535p_eta = create_diagram(FORMAT("decay N1535p to eta"), Decay_1_to_2, VMP,
                                                         {P.get("N1535p")},
                                                         {},
                                                         {Proton, P.get("eta")}
        );


        double const literature_value = P.get("N1535")->width() * P.get("N1535")->user_data<double>("branching_N_eta");

        Feynman_Process decay_Np({channel_decay_N1535p_eta});
        auto w1 = decay_Np.decay_width();
        auto const g = std::sqrt(literature_value / ( w1 ));
        couplings.set(coupl_str, g);
        buffer << FORMAT("{} {}\n", coupl_str, g);
    }
    {   /// N1535 to f0(500)
        auto coupl_str = coupling_string("N1535", "N", "f0_500");
        couplings.set(coupl_str, 1.);
        auto channel_decay_N1535p_f500_1 = create_diagram(FORMAT("decay N1535 to f_0(500) pi+pi-"), Decay_1_to_M2_1, VMP,
                                                          {P.get("N1535p")},
                                                          {P.get("f0_500")},
                                                          {Pi_Plus, Pi_Minus, Proton}
        );
        auto channel_decay_N1535p_f500_2 = create_diagram(FORMAT("decay N1535 to f_0(500) pi0pi0"), Decay_1_to_M2_1, VMP,
                                                          {P.get("N1535p")},
                                                          {P.get("f0_500")},
                                                          {Pi_Zero, Pi_Zero, Proton}
        );

        auto const literature_value = P.get("N1535")->width() * P.get("N1535")->user_data<double>("branching_N_pipi_f500");

        Feynman_Process decay_Np1({channel_decay_N1535p_f500_1});
        Feynman_Process decay_Np2({channel_decay_N1535p_f500_2});

        auto const w1 = decay_Np1.decay_width();
        auto const w2 = decay_Np2.decay_width();

        auto const g = std::sqrt(literature_value / ( w1 + w2 ));
        couplings.set(coupl_str, g);
        buffer << FORMAT("{} {}\n", coupl_str, g);
    }
    {   /// N1535 to D1232
        auto coupl_str = coupling_string("N1535", "D1232", "Pion");
        couplings.set(coupl_str, 1.);

        auto channel_decay_N1535p_D1232_pi_1 = create_diagram(FORMAT("decay N1535 to Delta pi 1"), Decay_1_to_M2_1, VMP,
                                                              {P.get("N1535p")},
                                                              {P.get("D1232pp")},
                                                              {Proton, Pi_Plus, Pi_Minus}
        );
        auto channel_decay_N1535p_D1232_pi_2 = create_diagram(FORMAT("decay N1535 to Delta pi 1"), Decay_1_to_M2_1_cross, VMP,
                                                              {P.get("N1535p")},
                                                              {P.get("D1232n")},
                                                              {Proton, Pi_Plus, Pi_Minus}
        );

        auto channel_decay_N1535p_D1232_pi_3 = create_diagram(FORMAT("decay N1535 to Delta pi 1"), Decay_1_to_M2_1, VMP,
                                                              {P.get("N1535p")},
                                                              {P.get("D1232p")},
                                                              {Proton, Pi_Zero, Pi_Zero}
        );
        auto channel_decay_N1535p_D1232_pi_4 = create_diagram(FORMAT("decay N1535 to Delta pi 1"), Decay_1_to_M2_1_cross, VMP,
                                                              {P.get("N1535p")},
                                                              {P.get("D1232p")},
                                                              {Proton, Pi_Zero, Pi_Zero}
        );

        auto channel_decay_N1535p_D1232_pi_5 = create_diagram(FORMAT("decay N1535 to Delta pi 1"), Decay_1_to_M2_1, VMP,
                                                              {P.get("N1535p")},
                                                              {P.get("D1232p")},
                                                              {Neutron, Pi_Plus, Pi_Zero}
        );
        auto channel_decay_N1535p_D1232_pi_6 = create_diagram(FORMAT("decay N1535 to Delta pi 1"), Decay_1_to_M2_1_cross, VMP,
                                                              {P.get("N1535p")},
                                                              {P.get("D1232n")},
                                                              {Neutron, Pi_Plus, Pi_Zero}
        );

        Feynman_Process decay_Np1({channel_decay_N1535p_D1232_pi_1, channel_decay_N1535p_D1232_pi_2});
        Feynman_Process decay_Np2({channel_decay_N1535p_D1232_pi_3, channel_decay_N1535p_D1232_pi_4});
        Feynman_Process decay_Np3({channel_decay_N1535p_D1232_pi_5, channel_decay_N1535p_D1232_pi_6});

        auto w1 = decay_Np1.decay_width();
        auto w2 = decay_Np2.decay_width();
        auto w3 = decay_Np3.decay_width();

        double const literature_value = P.get("N1535")->width() * P.get("N1535")->user_data<double>("branching_N_D1232_pi");

        auto const g = std::sqrt(literature_value / ( w1 + w2 + w3 ));
        couplings.set(coupl_str, g);
        buffer << FORMAT("{} {}\n", coupl_str, g);
    }
    {   /// N1535 to N1440
        auto coupl_str = coupling_string("N1535", "N1440", "Pion");
        couplings.set(coupl_str, 1.);

        auto channel_decay_N1535p_N1440_pi_2 = create_diagram(FORMAT("decay N1535 to Delta pi 1"), Decay_1_to_M2_1_cross, VMP,
                                                              {P.get("N1535p")},
                                                              {P.get("N1440n")},
                                                              {Proton, Pi_Plus, Pi_Minus}
        );

        auto channel_decay_N1535p_N1440_pi_3 = create_diagram(FORMAT("decay N1535 to Delta pi 1"), Decay_1_to_M2_1, VMP,
                                                              {P.get("N1535p")},
                                                              {P.get("N1440p")},
                                                              {Proton, Pi_Zero, Pi_Zero}
        );
        auto channel_decay_N1535p_N1440_pi_4 = create_diagram(FORMAT("decay N1535 to Delta pi 1"), Decay_1_to_M2_1_cross, VMP,
                                                              {P.get("N1535p")},
                                                              {P.get("N1440p")},
                                                              {Proton, Pi_Zero, Pi_Zero}
        );

        auto channel_decay_N1535p_N1440_pi_5 = create_diagram(FORMAT("decay N1535 to Delta pi 1"), Decay_1_to_M2_1, VMP,
                                                              {P.get("N1535p")},
                                                              {P.get("N1440p")},
                                                              {Neutron, Pi_Plus, Pi_Zero}
        );
        auto channel_decay_N1535p_N1440_pi_6 = create_diagram(FORMAT("decay N1535 to Delta pi 1"), Decay_1_to_M2_1_cross, VMP,
                                                              {P.get("N1535p")},
                                                              {P.get("N1440n")},
                                                              {Neutron, Pi_Plus, Pi_Zero}
        );

        auto channel_decay_N1535n_N1440_pi_5 = create_diagram(FORMAT("decay N1535 to Delta pi 1"), Decay_1_to_M2_1, VMP,
                                                              {P.get("N1535n")},
                                                              {P.get("N1440n")},
                                                              {Neutron, Pi_Zero, Pi_Zero}
        );
        auto channel_decay_N1535n_N1440_pi_6 = create_diagram(FORMAT("decay N1535 to Delta pi 1"), Decay_1_to_M2_1_cross, VMP,
                                                              {P.get("N1535n")},
                                                              {P.get("N1440n")},
                                                              {Neutron, Pi_Zero, Pi_Zero}
        );

        Feynman_Process decay_Np1({channel_decay_N1535p_N1440_pi_2});
        Feynman_Process decay_Np2({channel_decay_N1535p_N1440_pi_3, channel_decay_N1535p_N1440_pi_4});
        Feynman_Process decay_Np3({channel_decay_N1535p_N1440_pi_5, channel_decay_N1535p_N1440_pi_6});

        auto w1 = decay_Np1.decay_width();
        auto w2 = decay_Np2.decay_width();
        auto w3 = decay_Np3.decay_width();

        double const literature_value = P.get("N1535")->width() * P.get("N1535")->user_data<double>("branching_N_pipi_N1440");

        auto const g = std::sqrt(literature_value / ( w1 + w2 + w3 ));
        couplings.set(coupl_str, g);
        buffer << FORMAT("{} {}\n", coupl_str, g);
    }
    #endif
    #if D1600
    { // D1600 to N1440
        auto coupl_str = coupling_string("D1600", "N1440", "Pion");
        couplings.set(coupl_str, 1.);

        auto channel_decay_D1600_N1440_pin_pin_1 = create_diagram(FORMAT("decay D1600 to N1440 pi"), Decay_1_to_M2_1, VMP,
                                                              {P.get("D1600p")},
                                                              {P.get("N1440p")},
                                                              {Proton, Pi_Zero, Pi_Zero}
        );
        auto channel_decay_D1600_N1440_pin_pin_2 = create_diagram(FORMAT("decay D1600 to N1440 pi"), Decay_1_to_M2_1_cross, VMP,
                                                           {P.get("D1600p")},
                                                           {P.get("N1440p")},
                                                           {Proton, Pi_Zero, Pi_Zero}
        );

        auto channel_decay_D1600_N1440_pip_pin_1 = create_diagram(FORMAT("decay D1600 to N1440 pi"), Decay_1_to_M2_1, VMP,
                                                                  {P.get("D1600p")},
                                                                  {P.get("N1440p")},
                                                                  {Neutron, Pi_Plus, Pi_Zero}
        );
        auto channel_decay_D1600_N1440_pip_pin_2 = create_diagram(FORMAT("decay D1600 to N1440 pi"), Decay_1_to_M2_1_cross, VMP,
                                                                  {P.get("D1600p")},
                                                                  {P.get("N1440n")},
                                                                  {Neutron, Pi_Plus, Pi_Zero}
        );

        auto channel_decay_D1600_N1440_pip_pim_1 = create_diagram(FORMAT("decay D1600 to N1440 pi"), Decay_1_to_M2_1, VMP,
                                                                  {P.get("D1600p")},
                                                                  {P.get("N1440n")},
                                                                  {Proton, Pi_Minus, Pi_Plus}
        );

        Feynman_Process decay_Np1({channel_decay_D1600_N1440_pin_pin_1, channel_decay_D1600_N1440_pin_pin_2});
        Feynman_Process decay_Np2({channel_decay_D1600_N1440_pip_pin_1, channel_decay_D1600_N1440_pip_pin_2});
        Feynman_Process decay_Np3({channel_decay_D1600_N1440_pip_pim_1});
        auto w1 = decay_Np1.decay_width();
        auto w2 = decay_Np2.decay_width();
        auto w3 = decay_Np3.decay_width();
        double const literature_value = P.get("D1600")->width() * P.get("D1600")->user_data<double>("branching_N1440_pi");

        auto const g = std::sqrt(literature_value / ( w1 + w2 + w3 ));
        couplings.set(coupl_str, g);
        buffer << FORMAT("{} {}\n", coupl_str, g);
    };
    { // D1600 to D1232
        auto coupl_str = coupling_string("D1600", "D1232", "Pion");
        couplings.set(coupl_str, 1.);

        auto channel_decay_D1600_D1232_pip_pin_1 = create_diagram(FORMAT("decay D1600 to D1232 pi+ pi0 1"), Decay_1_to_M2_1, VMP,
                                                           {P.get("D1600pp")},
                                                           {P.get("D1232pp")},
                                                           {Proton, Pi_Plus, Pi_Zero}
        );
        auto channel_decay_D1600_D1232_pip_pin_2 = create_diagram(FORMAT("decay D1600 to D1232 pi+ pi0 2"), Decay_1_to_M2_1_cross, VMP,
                                                           {P.get("D1600pp")},
                                                           {P.get("D1232p")},
                                                           {Proton, Pi_Plus, Pi_Zero}
        );

        auto channel_decay_D1600_D1232_pip_pip_1 = create_diagram(FORMAT("decay D1600 to D1232 pi+ pi+ 1"), Decay_1_to_M2_1, VMP,
                                                                  {P.get("D1600pp")},
                                                                  {P.get("D1232p")},
                                                                  {Neutron, Pi_Plus, Pi_Plus}
        );
        auto channel_decay_D1600_D1232_pip_pip_2 = create_diagram(FORMAT("decay D1600 to D1232 pi+ pi+ 2"), Decay_1_to_M2_1_cross, VMP,
                                                                  {P.get("D1600pp")},
                                                                  {P.get("D1232p")},
                                                                  {Neutron, Pi_Plus, Pi_Plus}
        );

        Feynman_Process decay_Np1({channel_decay_D1600_D1232_pip_pin_1, channel_decay_D1600_D1232_pip_pin_2});
        Feynman_Process decay_Np2({channel_decay_D1600_D1232_pip_pip_1, channel_decay_D1600_D1232_pip_pip_2});
        auto w1 = decay_Np1.decay_width();
        auto w2 = decay_Np2.decay_width();

        double const literature_value = P.get("D1600")->width() * P.get("D1600")->user_data<double>("branching_D1232_pi_pi");

        auto const g = std::sqrt(literature_value / ( w1 + w2 ));
        couplings.set(coupl_str, g);
        buffer << FORMAT("{} {}\n", coupl_str, g);
    }
    #endif
    #if D1620
    { // D1620 to D1232
        auto coupl_str = coupling_string("D1620", "D1232", "Pion");
        couplings.set(coupl_str, 1.);

        auto channel_decay_D1620_D1232_pip_pin_1 = create_diagram(FORMAT("decay D1620 to D1232 pi+ pi0 1"), Decay_1_to_M2_1, VMP,
                                                                  {P.get("D1620pp")},
                                                                  {P.get("D1232pp")},
                                                                  {Proton, Pi_Plus, Pi_Zero}
        );
        auto channel_decay_D1620_D1232_pip_pin_2 = create_diagram(FORMAT("decay D1620 to D1232 pi+ pi0 2"), Decay_1_to_M2_1_cross, VMP,
                                                                  {P.get("D1620pp")},
                                                                  {P.get("D1232p")},
                                                                  {Proton, Pi_Plus, Pi_Zero}
        );

        auto channel_decay_D1620_D1232_pip_pip_1 = create_diagram(FORMAT("decay D1620 to D1232 pi+ pi+ 1"), Decay_1_to_M2_1, VMP,
                                                                  {P.get("D1620pp")},
                                                                  {P.get("D1232p")},
                                                                  {Neutron, Pi_Plus, Pi_Plus}
        );
        auto channel_decay_D1620_D1232_pip_pip_2 = create_diagram(FORMAT("decay D1620 to D1232 pi+ pi+ 2"), Decay_1_to_M2_1_cross, VMP,
                                                                  {P.get("D1620pp")},
                                                                  {P.get("D1232p")},
                                                                  {Neutron, Pi_Plus, Pi_Plus}
        );

        Feynman_Process decay_Np1({channel_decay_D1620_D1232_pip_pin_1, channel_decay_D1620_D1232_pip_pin_2});
        Feynman_Process decay_Np2({channel_decay_D1620_D1232_pip_pip_1, channel_decay_D1620_D1232_pip_pip_2});
        auto w1 = decay_Np1.decay_width();
        auto w2 = decay_Np2.decay_width();

        double const literature_value = P.get("D1620")->width() * P.get("D1620")->user_data<double>("branching_D1232_pi_pi");

        auto const g = std::sqrt(literature_value / ( w1 + w2 ));
        couplings.set(coupl_str, g);
        buffer << FORMAT("{} {}\n", coupl_str, g);
    }
    { // D1620 to N1440
        auto coupl_str = coupling_string("D1620", "N1440", "Pion");
        couplings.set(coupl_str, 1.);

        auto channel_decay_D1620_N1440_pin_pin_1 = create_diagram(FORMAT("decay D1620 to N1440 pi"), Decay_1_to_M2_1, VMP,
                                                                  {P.get("D1620p")},
                                                                  {P.get("N1440p")},
                                                                  {Proton, Pi_Zero, Pi_Zero}
        );
        auto channel_decay_D1620_N1440_pin_pin_2 = create_diagram(FORMAT("decay D1620 to N1440 pi"), Decay_1_to_M2_1_cross, VMP,
                                                                  {P.get("D1620p")},
                                                                  {P.get("N1440p")},
                                                                  {Proton, Pi_Zero, Pi_Zero}
        );

        auto channel_decay_D1620_N1440_pip_pin_1 = create_diagram(FORMAT("decay D1620 to N1440 pi"), Decay_1_to_M2_1, VMP,
                                                                  {P.get("D1620p")},
                                                                  {P.get("N1440p")},
                                                                  {Neutron, Pi_Plus, Pi_Zero}
        );
        auto channel_decay_D1620_N1440_pip_pin_2 = create_diagram(FORMAT("decay D1620 to N1440 pi"), Decay_1_to_M2_1_cross, VMP,
                                                                  {P.get("D1620p")},
                                                                  {P.get("N1440n")},
                                                                  {Neutron, Pi_Plus, Pi_Zero}
        );

        auto channel_decay_D1620_N1440_pip_pim_1 = create_diagram(FORMAT("decay D1620 to N1440 pi"), Decay_1_to_M2_1, VMP,
                                                                  {P.get("D1620p")},
                                                                  {P.get("N1440n")},
                                                                  {Proton, Pi_Minus, Pi_Plus}
        );

        Feynman_Process decay_Np1({channel_decay_D1620_N1440_pin_pin_1, channel_decay_D1620_N1440_pin_pin_2});
        Feynman_Process decay_Np2({channel_decay_D1620_N1440_pip_pin_1, channel_decay_D1620_N1440_pip_pin_2});
        Feynman_Process decay_Np3({channel_decay_D1620_N1440_pip_pim_1});
        auto w1 = decay_Np1.decay_width();
        auto w2 = decay_Np2.decay_width();
        auto w3 = decay_Np3.decay_width();
        double const literature_value = P.get("D1620")->width() * P.get("D1620")->user_data<double>("branching_N1440_pi");

        auto const g = std::sqrt(literature_value / ( w1 + w2 + w3 ));
        couplings.set(coupl_str, g);
        buffer << FORMAT("{} {}\n", coupl_str, g);
    };
    { // D1620 to Rho
        auto coupl_str = coupling_string("D1620", "Rho", "N");
        couplings.set(coupl_str, 1.);

        auto channel_decay_D1620_rho_pip_pin_1 = create_diagram(FORMAT("decay D1620 to D1232 pi+ pi0 1"), Decay_1_to_1_M2, VMP,
                                                                  {P.get("D1620pp")},
                                                                  {P.get("rho+")},
                                                                  {Proton, Pi_Plus, Pi_Zero}
        );

        Feynman_Process decay_Np1({channel_decay_D1620_rho_pip_pin_1});
        auto w1 = decay_Np1.decay_width();

        double const literature_value = P.get("D1620")->width() * P.get("D1620")->user_data<double>("branching_rho_pi_pi");

        auto const g = std::sqrt(literature_value / ( w1 ));
        couplings.set(coupl_str, g);
        buffer << FORMAT("{} {}\n", coupl_str, g);
    }
    #endif
    #if D1700
    { // D1700 to D1232
        auto coupl_str = coupling_string("D1700", "D1232", "Pion");
        couplings.set(coupl_str, 1.);

        auto channel_decay_D1700_D1232_pip_pin_1 = create_diagram(FORMAT("decay D1700 to D1232 pi+ pi0 1"), Decay_1_to_M2_1, VMP,
                                                                  {P.get("D1700pp")},
                                                                  {P.get("D1232pp")},
                                                                  {Proton, Pi_Plus, Pi_Zero}
        );
        auto channel_decay_D1700_D1232_pip_pin_2 = create_diagram(FORMAT("decay D1700 to D1232 pi+ pi0 2"), Decay_1_to_M2_1_cross, VMP,
                                                                  {P.get("D1700pp")},
                                                                  {P.get("D1232p")},
                                                                  {Proton, Pi_Plus, Pi_Zero}
        );

        auto channel_decay_D1700_D1232_pip_pip_1 = create_diagram(FORMAT("decay D1700 to D1232 pi+ pi+ 1"), Decay_1_to_M2_1, VMP,
                                                                  {P.get("D1700pp")},
                                                                  {P.get("D1232p")},
                                                                  {Neutron, Pi_Plus, Pi_Plus}
        );
        auto channel_decay_D1700_D1232_pip_pip_2 = create_diagram(FORMAT("decay D1700 to D1232 pi+ pi+ 2"), Decay_1_to_M2_1_cross, VMP,
                                                                  {P.get("D1700pp")},
                                                                  {P.get("D1232p")},
                                                                  {Neutron, Pi_Plus, Pi_Plus}
        );

        Feynman_Process decay_Np1({channel_decay_D1700_D1232_pip_pin_1, channel_decay_D1700_D1232_pip_pin_2});
        Feynman_Process decay_Np2({channel_decay_D1700_D1232_pip_pip_1, channel_decay_D1700_D1232_pip_pip_2});
        auto w1 = decay_Np1.decay_width();
        auto w2 = decay_Np2.decay_width();

        double const literature_value = P.get("D1700")->width() * P.get("D1700")->user_data<double>("branching_D1232_pi_pi");

        auto const g = std::sqrt(literature_value / ( w1 + w2 ));
        couplings.set(coupl_str, g);
        buffer << FORMAT("{} {}\n", coupl_str, g);
    }
    { // D1700 to N1440
        auto coupl_str = coupling_string("D1700", "N1440", "Pion");
        couplings.set(coupl_str, 1.);

        auto channel_decay_D1700_N1440_pin_pin_1 = create_diagram(FORMAT("decay D1700 to N1440 pi"), Decay_1_to_M2_1, VMP,
                                                                  {P.get("D1700p")},
                                                                  {P.get("N1440p")},
                                                                  {Proton, Pi_Zero, Pi_Zero}
        );
        auto channel_decay_D1700_N1440_pin_pin_2 = create_diagram(FORMAT("decay D1700 to N1440 pi"), Decay_1_to_M2_1_cross, VMP,
                                                                  {P.get("D1700p")},
                                                                  {P.get("N1440p")},
                                                                  {Proton, Pi_Zero, Pi_Zero}
        );

        auto channel_decay_D1700_N1440_pip_pin_1 = create_diagram(FORMAT("decay D1700 to N1440 pi"), Decay_1_to_M2_1, VMP,
                                                                  {P.get("D1700p")},
                                                                  {P.get("N1440p")},
                                                                  {Neutron, Pi_Plus, Pi_Zero}
        );
        auto channel_decay_D1700_N1440_pip_pin_2 = create_diagram(FORMAT("decay D1700 to N1440 pi"), Decay_1_to_M2_1_cross, VMP,
                                                                  {P.get("D1700p")},
                                                                  {P.get("N1440n")},
                                                                  {Neutron, Pi_Plus, Pi_Zero}
        );

        auto channel_decay_D1700_N1440_pip_pim_1 = create_diagram(FORMAT("decay D1700 to N1440 pi"), Decay_1_to_M2_1, VMP,
                                                                  {P.get("D1700p")},
                                                                  {P.get("N1440n")},
                                                                  {Proton, Pi_Minus, Pi_Plus}
        );

        Feynman_Process decay_Np1({channel_decay_D1700_N1440_pin_pin_1, channel_decay_D1700_N1440_pin_pin_2});
        Feynman_Process decay_Np2({channel_decay_D1700_N1440_pip_pin_1, channel_decay_D1700_N1440_pip_pin_2});
        Feynman_Process decay_Np3({channel_decay_D1700_N1440_pip_pim_1});
        auto w1 = decay_Np1.decay_width();
        auto w2 = decay_Np2.decay_width();
        auto w3 = decay_Np3.decay_width();
        double const literature_value = P.get("D1700")->width() * P.get("D1700")->user_data<double>("branching_N1440_pi");

        auto const g = std::sqrt(literature_value / ( w1 + w2 + w3 ));
        couplings.set(coupl_str, g);
        buffer << FORMAT("{} {}\n", coupl_str, g);
    };
    { // D1700 to Rho
        auto coupl_str = coupling_string("D1700", "Rho", "N");
        couplings.set(coupl_str, 1.);

        auto channel_decay_D1700_rho_pip_pin_1 = create_diagram(FORMAT("decay D1700 to D1232 pi+ pi0 1"), Decay_1_to_1_M2, VMP,
                                                                {P.get("D1700pp")},
                                                                {P.get("rho+")},
                                                                {Proton, Pi_Plus, Pi_Zero}
        );

        Feynman_Process decay_Np1({channel_decay_D1700_rho_pip_pin_1});
        auto w1 = decay_Np1.decay_width();

        double const literature_value = P.get("D1700")->width() * P.get("D1700")->user_data<double>("branching_rho_pi_pi");

        auto const g = std::sqrt(literature_value / ( w1 ));
        couplings.set(coupl_str, g);
        buffer << FORMAT("{} {}\n", coupl_str, g);
    }
    #endif
    std::ofstream coupling_constants_out("./data/coupling_constants_isospin_symmetry.txt");
    coupling_constants_out << buffer.str();
	return EXIT_SUCCESS;
}