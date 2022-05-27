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
	auto const& Proton = P.get("proton");
	auto const& Neutron = P.get("neutron");
	//auto const& Photon = P.get("photon");
	auto const& Pi_Zero = P.get("pi0");
	auto const& Pi_Minus = P.get("pi-");
	auto const& Pi_Plus = P.get("pi+");

	for( auto& [key, particle] : P ){
		auto k = key;
		auto p = particle;
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

	/* Other */
    {   /// f0(500) to pi pi
        auto particle = P.get("f0_500");
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
        std::cout << FORMAT("{} -> {} {} g: {}\n", particle->name(), Pi_Plus->name(), Pi_Minus->name(), g);
    }
    {   /// N1440 to D1232
        P.get("N1440p")->user_data("form_factor", identity);
        P.get("N1440p")->user_data("gD1232", 1.0);
        P.get("N1440n")->user_data("form_factor", identity);
	    P.get("N1440n")->user_data("gD1232", 1.0);
        P.get("D1232pp")->user_data("form_factor", identity);
        P.get("D1232p")->user_data("form_factor", identity);
        P.get("D1232n")->user_data("form_factor", identity);
        P.get("D1232m")->user_data("form_factor", identity);

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

        auto channel_decay_N1440n_D1232_pi_1 = create_diagram(FORMAT("decay N1440 to Delta pi 1"), Decay_1_to_M2_1, VMP,
                                                              {P.get("N1440n")},
                                                              {P.get("D1232n")},
                                                              {Proton, Pi_Minus, Pi_Zero}
        );
        auto channel_decay_N1440n_D1232_pi_2 = create_diagram(FORMAT("decay N1440 to Delta pi 1"), Decay_1_to_M2_1_cross, VMP,
                                                              {P.get("N1440n")},
                                                              {P.get("D1232p")},
                                                              {Proton, Pi_Minus, Pi_Zero}
        );

        auto channel_decay_N1440n_D1232_pi_3 = create_diagram(FORMAT("decay N1440 to Delta pi 1"), Decay_1_to_M2_1, VMP,
                                                         {P.get("N1440n")},
                                                         {P.get("D1232p")},
                                                         {Neutron, Pi_Plus, Pi_Minus}
        );
        auto channel_decay_N1440n_D1232_pi_4 = create_diagram(FORMAT("decay N1440 to Delta pi 1"), Decay_1_to_M2_1_cross, VMP,
                                                         {P.get("N1440n")},
                                                         {P.get("D1232m")},
                                                         {Neutron, Pi_Plus, Pi_Minus}
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

        Feynman_Process decay_Nn1({channel_decay_N1440n_D1232_pi_1, channel_decay_N1440n_D1232_pi_2});
        Feynman_Process decay_Nn2({channel_decay_N1440n_D1232_pi_3, channel_decay_N1440n_D1232_pi_4});
        Feynman_Process decay_Nn3({channel_decay_N1440n_D1232_pi_5, channel_decay_N1440n_D1232_pi_6});

        auto w1 = decay_Np1.decay_width();
        auto w2 = decay_Np2.decay_width();
        auto w3 = decay_Np3.decay_width();
        auto w4 = decay_Nn1.decay_width();
        auto w5 = decay_Nn2.decay_width();
        auto w6 = decay_Nn3.decay_width();

        double const literature_value = P.get("N1440")->width() * ( P.get("N1440")->user_data<double>("branching_N_pipi_D1232_upper") + P.get("N1440")->user_data<double>("branching_N_pipi_D1232_lower")) / 2.;

        std::cout << FORMAT("g(N1440+ -> D1232): " ) << std::setw(10) << std::setprecision(10)<< std::sqrt(literature_value / ( w1 + w2 + w3 )) << "\n";
        std::cout << FORMAT("g(N1440n -> D1232): " ) << std::setw(10) << std::setprecision(10)<< std::sqrt(literature_value / ( w4 + w5 + w6 )) << "\n";
    }
	{   /// N1440 to f0(500)
		P.get("N1440p")->user_data("gRNf0_500", 1.0);
		P.get("N1440n")->user_data("gRNf0_500", 1.0);
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

		auto channel_decay_N1440n_f500_1 = create_diagram(FORMAT("decay N1440 to f_0(500) pi+pi-"), Decay_1_to_M2_1, VMP,
		                                                  {P.get("N1440n")},
		                                                  {P.get("f0_500")},
		                                                  {Pi_Plus, Pi_Minus, Neutron}
		);
		auto channel_decay_N1440n_f500_2 = create_diagram(FORMAT("decay N1440 to f_0(500) pi0pi0"), Decay_1_to_M2_1, VMP,
		                                                  {P.get("N1440n")},
		                                                  {P.get("f0_500")},
		                                                  {Pi_Zero, Pi_Zero, Neutron}
		);

		double const literature_value = P.get("N1440")->width() * ( P.get("N1440")->user_data<double>("branching_N_pipi_D1232_upper") + P.get("N1440")->user_data<double>("branching_N_pipi_D1232_lower")) / 2.;

		Feynman_Process decay_Np1({channel_decay_N1440p_f500_1});
		Feynman_Process decay_Np2({channel_decay_N1440p_f500_2});
		Feynman_Process decay_Nn1({channel_decay_N1440n_f500_1});
		Feynman_Process decay_Nn2({channel_decay_N1440n_f500_2});
		auto w1 = decay_Np1.decay_width();
		auto w2 = decay_Np2.decay_width();
		auto w3 = decay_Nn1.decay_width();
		auto w4 = decay_Nn2.decay_width();

		std::cout << FORMAT("g(N1440+ -> f_0(500): " ) << std::setw(10) << std::setprecision(10)<< std::sqrt(literature_value / ( w1 + w2 )) << "\n";
		std::cout << FORMAT("g(N1440n -> f_0(500): " ) << std::setw(10) << std::setprecision(10)<< std::sqrt(literature_value / ( w3 + w4 )) << "\n";
	}
	{   /// N1520 to rho
		P.get("N1520p")->user_data("gRNrho", 1.0);
		P.get("N1520n")->user_data("gRNrho", 1.0);
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

		auto channel_decay_N1520n_rho_1 = create_diagram(FORMAT("decay N1520n to rho0 pi+pi-"), Decay_1_to_M2_1, VMP,
		                                                 {P.get("N1520n")},
		                                                 {P.get("rho0")},
		                                                 {Pi_Plus, Pi_Minus, Neutron}
		);
		auto channel_decay_N1520n_rho_2 = create_diagram(FORMAT("decay N1520n to rho0 pi-pi0"), Decay_1_to_M2_1, VMP,
		                                                 {P.get("N1520n")},
		                                                 {P.get("rho-")},
		                                                 {Pi_Minus, Pi_Zero, Proton}
		);

		double const literature_value = P.get("N1520")->width() * ( P.get("N1520")->user_data<double>("branching_N_pipi_rho_upper") + P.get("N1520")->user_data<double>("branching_N_pipi_rho_lower")) / 2.;

		Feynman_Process decay_Np1({channel_decay_N1520p_rho_1});
		Feynman_Process decay_Np2({channel_decay_N1520p_rho_2});
		Feynman_Process decay_Nn1({channel_decay_N1520n_rho_1});
		Feynman_Process decay_Nn2({channel_decay_N1520n_rho_2});
		auto w1 = decay_Np1.decay_width();
		auto w2 = decay_Np2.decay_width();
		auto w4 = decay_Nn1.decay_width();
		auto w5 = decay_Nn2.decay_width();

		std::cout << FORMAT("g(N1520+ -> rho: " ) << std::setw(10) << std::setprecision(10)<< std::sqrt(literature_value / ( w1 + w2  )) << "\n";
		std::cout << FORMAT("g(N1520n -> rho: " ) << std::setw(10) << std::setprecision(10)<< std::sqrt(literature_value / ( w4 + w5  )) << "\n";
	}
	{   /// N1520 to D1232
		P.get("N1520p")->user_data("form_factor", identity);
		P.get("N1520p")->user_data("gD1232", 1.0);
		P.get("N1520n")->user_data("form_factor", identity);
		P.get("N1520n")->user_data("gD1232", 1.0);
		P.get("D1232pp")->user_data("form_factor", identity);
		P.get("D1232p")->user_data("form_factor", identity);
		P.get("D1232n")->user_data("form_factor", identity);
		P.get("D1232m")->user_data("form_factor", identity);

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

		auto channel_decay_N1520n_D1232_pi_1 = create_diagram(FORMAT("decay N1520 to Delta pi 1"), Decay_1_to_M2_1, VMP,
		                                                      {P.get("N1520n")},
		                                                      {P.get("D1232n")},
		                                                      {Proton, Pi_Minus, Pi_Zero}
		);
		auto channel_decay_N1520n_D1232_pi_2 = create_diagram(FORMAT("decay N1520 to Delta pi 1"), Decay_1_to_M2_1_cross, VMP,
		                                                      {P.get("N1520n")},
		                                                      {P.get("D1232p")},
		                                                      {Proton, Pi_Minus, Pi_Zero}
		);

		auto channel_decay_N1520n_D1232_pi_3 = create_diagram(FORMAT("decay N1520 to Delta pi 1"), Decay_1_to_M2_1, VMP,
		                                                      {P.get("N1520n")},
		                                                      {P.get("D1232p")},
		                                                      {Neutron, Pi_Plus, Pi_Minus}
		);
		auto channel_decay_N1520n_D1232_pi_4 = create_diagram(FORMAT("decay N1520 to Delta pi 1"), Decay_1_to_M2_1_cross, VMP,
		                                                      {P.get("N1520n")},
		                                                      {P.get("D1232m")},
		                                                      {Neutron, Pi_Plus, Pi_Minus}
		);

		auto channel_decay_N1520n_D1232_pi_5 = create_diagram(FORMAT("decay N1520 to Delta pi 1"), Decay_1_to_M2_1, VMP,
		                                                      {P.get("N1520n")},
		                                                      {P.get("D1232n")},
		                                                      {Neutron, Pi_Zero, Pi_Zero}
		);
		auto channel_decay_N1520n_D1232_pi_6 = create_diagram(FORMAT("decay N1520 to Delta pi 1"), Decay_1_to_M2_1_cross, VMP,
		                                                      {P.get("N1520n")},
		                                                      {P.get("D1232n")},
		                                                      {Neutron, Pi_Zero, Pi_Zero}
		);

		Feynman_Process decay_Np1({channel_decay_N1520p_D1232_pi_1, channel_decay_N1520p_D1232_pi_2});
		Feynman_Process decay_Np2({channel_decay_N1520p_D1232_pi_3, channel_decay_N1520p_D1232_pi_4});
		Feynman_Process decay_Np3({channel_decay_N1520p_D1232_pi_5, channel_decay_N1520p_D1232_pi_6});

		Feynman_Process decay_Nn1({channel_decay_N1520n_D1232_pi_1, channel_decay_N1520n_D1232_pi_2});
		Feynman_Process decay_Nn2({channel_decay_N1520n_D1232_pi_3, channel_decay_N1520n_D1232_pi_4});
		Feynman_Process decay_Nn3({channel_decay_N1520n_D1232_pi_5, channel_decay_N1520n_D1232_pi_6});

		auto w1 = decay_Np1.decay_width();
		auto w2 = decay_Np2.decay_width();
		auto w3 = decay_Np3.decay_width();
		auto w4 = decay_Nn1.decay_width();
		auto w5 = decay_Nn2.decay_width();
		auto w6 = decay_Nn3.decay_width();

		double const literature_value = P.get("N1520")->width() * ( P.get("N1520")->user_data<double>("branching_N_pipi_D1232_upper") + P.get("N1520")->user_data<double>("branching_N_pipi_D1232_lower")) / 2.;

		std::cout << FORMAT("g(N1520+ -> D1232): " ) << std::setw(10) << std::setprecision(10)<< std::sqrt(literature_value / ( w1 + w2 + w3 )) << "\n";
		std::cout << FORMAT("g(N1520n -> D1232): " ) << std::setw(10) << std::setprecision(10)<< std::sqrt(literature_value / ( w4 + w5 + w6 )) << "\n";
	}

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
		auto& Dm = particles_Dm[i];
		if( Dpp->spin().j() > 2 ){
			continue;
		}
		Dpp->user_data("gRNpi", 1.);
		Dpp->user_data("form_factor", identity);
		Dm->user_data("gRNpi", 1.);
		Dm->user_data("form_factor", identity);

		auto channel_decay_Dpp1 = create_diagram(FORMAT("decay {} to proton pi+", Dpp->name()), Decay_1_to_2, VMP,
		                                             {Dpp},
		                                             {},
		                                             {Proton, Pi_Plus}
		);

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
		std::cout << FORMAT("g({} -> N + pi): ", Dm->name())  << std::setw(10) << std::setprecision(10) << std::sqrt(literature_value4 / ( w6 )) << "\n";
	}
	return EXIT_SUCCESS;
}