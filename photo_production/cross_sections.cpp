#include <feynumeric/command_line_manager.hpp>
#include <feynumeric/core.hpp>
#include <feynumeric/qed.hpp>


#include "effective_lagrangian_model.hpp"
#include "form_factors.hpp"

#include <iostream>

int main(int argc, char** argv)
{
	using namespace Feynumeric;
	using namespace Feynumeric::Units;

	std::string const CMD_PROCESS_ELASTIC_SCATTERING = "elasticscattering";
	std::string const CMD_PROCESS_PHOTO_PRODUCTION   = "photoproduction";

	Command_Line_Manager cmd(argc, argv);

	cmd.register_command("particle_file", true, "file with particle parameters");
	cmd.register_command("form_factor", std::string("none"), "which form factor to use (none, cassing, manley, moniz, cutkosky)");
	cmd.register_command("channel", std::string("s"), "which channel to use [s, t, u, c] or any combination.");
	cmd.register_command("process", std::string("photoproduction"), FORMAT("which process to use: {} {}", CMD_PROCESS_PHOTO_PRODUCTION, CMD_PROCESS_ELASTIC_SCATTERING));
	cmd.register_command("sqrt_s", true, "the energy in the center of mass frame in GeV");
	cmd.register_command("help", false, "list all command line parameters");

	if( cmd.is_enabled("help") ){
		return EXIT_SUCCESS;
	}

	cmd.crash_on_missing_mandatory_command();

	double sqrt_s = cmd.as_double("sqrt_s");

	auto const& channel = cmd.as_string("channel");
	bool const s_channel = channel.find('s') != std::string::npos;
	bool const u_channel = channel.find('u') != std::string::npos;
	bool const t_channel = channel.find('t') != std::string::npos;
	bool const c_channel = channel.find('c') != std::string::npos;




	Particle_Manager P(cmd.as_string("particle_file"));
	Particle_Ptr const& D1232pp  = P["D1232pp"];
	Particle_Ptr const& D1232p   = P["D1232p"];
	Particle_Ptr const& D1232n   = P["D1232n"];
	Particle_Ptr const& D1232m   = P["D1232m"];
	Particle_Ptr const& N1440p   = P["N1440p"];
	Particle_Ptr const& N1440n   = P["N1440n"];
	Particle_Ptr const& N1520p   = P["N1520p"];
	Particle_Ptr const& N1520n   = P["N1520n"];
	Particle_Ptr const& Proton   = P["proton"];
	Particle_Ptr const& Neutron  = P["neutron"];
	Particle_Ptr const& Pi_Plus  = P["pi+"];
	Particle_Ptr const& Pi_Minus = P["pi-"];
	Particle_Ptr const& Pi_Zero  = P["pi0"];

	std::vector<Particle_Ptr> resonances = {D1232pp, D1232p, D1232n, D1232m, N1440p, N1440n, N1520p, N1520n};

	std::string const& form_factor = cmd.as_string("form_factor");

	FORM_FACTOR_FUNCTION ff;

	if( form_factor == "none" ){
		ff = identity;
	}
	else if( form_factor == "moniz" ){
		ff = moniz;
	}
	else if( form_factor == "manley" ){
		ff = manley;
	}
	else if( form_factor == "cutkostky" ){
		ff = cutkosky;
	}
	else if( form_factor == "cassing" ){
		ff = cassing;
	}
	else{
		critical_error("Unknown form factor");
	}

	for( auto& resonance : resonances ){
		resonance->user_data("form_factor", ff);
		resonance->width([&](double p2){
			double const beta = 300._MeV;
			double const q0 = momentum(resonance->mass(), Neutron->mass(), Pi_Plus->mass());
			double const q = momentum(std::sqrt(p2), Neutron->mass(), Pi_Plus->mass());
			int l = resonance->user_data<double>("l");
			return resonance->width() * ff(beta, q0, q, l);
		});
	}

	if( !cmd.is_enabled("D1232") &&  !cmd.is_enabled("N1440") && !cmd.is_enabled("N1520") && !cmd.is_enabled("N1535") ){
		critical_error("No particle activated.\n");
	}

	init_vertices(P);

	if( cmd.as_string("process") == CMD_PROCESS_ELASTIC_SCATTERING ){
		std::vector<Feynman_Diagram_Ptr> diagrams_proton_pi_plus;
		if( cmd.is_enabled("D1232") ){
			if( s_channel ){
				diagrams_proton_pi_plus.push_back(
						create_diagram("D1232 s", Double_Wrench, VMP,
						               {Proton, Pi_Plus},
						               {D1232pp},
						               {Pi_Plus, Proton}
						));
			}
		}
		Feynman_Process scattering_proton_pi_plus(diagrams_proton_pi_plus);

		status("Scattering proton pi_plus -> proton pi_plus");
		status(FORMAT("sqrt_s = {}", sqrt_s));
		scattering_proton_pi_plus.print_dsigma_dcos_table(std::cout, sqrt_s, 0.1);
		scattering_proton_pi_plus.print_sigma_table(std::cout, {sqrt_s});
		std::cout << "\n\n";
	}
	else if( cmd.as_string("process") == CMD_PROCESS_PHOTO_PRODUCTION ){
		std::vector<Feynman_Diagram_Ptr> diagrams_proton_pi_zero;
		std::vector<Feynman_Diagram_Ptr> diagrams_proton_pi_minus;
		std::vector<Feynman_Diagram_Ptr> diagrams_neutron_pi_plus;



		if( cmd.is_enabled("N1440") )
		{
			if( s_channel ){
				diagrams_proton_pi_zero.push_back(
						create_diagram("N1440 s", Double_Wrench, VMP,
						               {Proton, QED::Photon},
						               {N1440p},
						               {Pi_Zero, Proton}
						));
				diagrams_proton_pi_minus.push_back(
						create_diagram("N1440 s", Double_Wrench, VMP,
						               {Neutron, QED::Photon},
						               {N1440n},
						               {Pi_Minus, Proton}
						));
				diagrams_neutron_pi_plus.push_back(
						create_diagram("N1440 s", Double_Wrench, VMP,
						               {Proton, QED::Photon},
						               {N1440p},
						               {Pi_Plus, Neutron}
						));
			}
			if( u_channel ){
				diagrams_proton_pi_zero.push_back(
						create_diagram("N1440 u", X_Man, VMP,
						               {Proton, QED::Photon},
						               {N1440p},
						               {Pi_Zero, Proton}
						));
				diagrams_proton_pi_minus.push_back(
						create_diagram("N1440 u", X_Man, VMP,
						               {Neutron, QED::Photon},
						               {N1440p},
						               {Pi_Minus, Proton}
						));
				diagrams_neutron_pi_plus.push_back(
						create_diagram("N1440 u", X_Man, VMP,
						               {Proton, QED::Photon},
						               {N1440n},
						               {Pi_Plus, Neutron}
						));
			}
		}
		if( cmd.is_enabled("N1520") )
		{
			if( s_channel ){
				diagrams_proton_pi_zero.push_back(
						create_diagram("N1520 s", Double_Wrench, VMP,
						               {Proton, QED::Photon},
						               {N1520p},
						               {Pi_Zero, Proton}
						));
				diagrams_proton_pi_minus.push_back(
						create_diagram("N1520 s", Double_Wrench, VMP,
						               {Neutron, QED::Photon},
						               {N1520n},
						               {Pi_Minus, Proton}
						));
				diagrams_neutron_pi_plus.push_back(
						create_diagram("N1520 s", Double_Wrench, VMP,
						               {Proton, QED::Photon},
						               {N1520p},
						               {Pi_Plus, Neutron}
						));
			}
			if( u_channel ){
				diagrams_proton_pi_zero.push_back(
						create_diagram("N1520 u", X_Man, VMP,
						               {Proton, QED::Photon},
						               {N1520p},
						               {Pi_Zero, Proton}
						));
				diagrams_proton_pi_minus.push_back(
						create_diagram("N1520 u", X_Man, VMP,
						               {Neutron, QED::Photon},
						               {N1520p},
						               {Pi_Minus, Proton}
						));
				diagrams_neutron_pi_plus.push_back(
						create_diagram("N1520 u", X_Man, VMP,
						               {Proton, QED::Photon},
						               {N1520n},
						               {Pi_Plus, Neutron}
						));
			}
		}

		Feynman_Process scattering_proton_pi_zero(diagrams_proton_pi_zero);
		Feynman_Process scattering_proton_pi_minus(diagrams_proton_pi_minus);
		Feynman_Process scattering_neutron_pi_plus(diagrams_neutron_pi_plus);

		status("Scattering proton gamma -> proton pi0");
		status(FORMAT("sqrt_s = {}", sqrt_s));
		scattering_proton_pi_zero.print_dsigma_dcos_table(std::cout, sqrt_s, 0.1);
	//	scattering_proton_pi_zero.print_sigma_table(std::cout, {{1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0}});
		std::cout << "\n\n";
		status("Scattering neutron gamma -> proton pi-");
		status(FORMAT("sqrt_s = {}", sqrt_s));
		scattering_proton_pi_minus.print_dsigma_dcos_table(std::cout, sqrt_s, 0.1);
	//	scattering_proton_pi_minus.print_sigma_table(std::cout, {{1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0}});
		status("Scattering proton gamma -> neutron pi+");
		status(FORMAT("sqrt_s = {}", sqrt_s));
		scattering_neutron_pi_plus.print_dsigma_dcos_table(std::cout, sqrt_s, 0.1);
	//	scattering_neutron_pi_plus.print_sigma_table(std::cout, {{1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0}});
	}
	return EXIT_SUCCESS;
}