#include <feynumeric/command_line_manager.hpp>
#include <feynumeric/core.hpp>
#include <feynumeric/qed.hpp>
#include <feynumeric/timer.hpp>

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
	cmd.register_command("sqrt_s", false, "the energy in the center of mass frame in GeV");
	cmd.register_command("help", false, "list all command line parameters");

	if( cmd.is_enabled("help") ){
		return EXIT_SUCCESS;
	}

	cmd.crash_on_missing_mandatory_command();

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
	Particle_Ptr const& D1600pp  = P["D1600pp"];
	Particle_Ptr const& D1600p   = P["D1600p"];
	Particle_Ptr const& D1600n   = P["D1600n"];
	Particle_Ptr const& D1600m   = P["D1600m"];
	Particle_Ptr const& D1920pp  = P["D1920pp"];
	Particle_Ptr const& D1920p   = P["D1920p"];
	Particle_Ptr const& D1920n   = P["D1920n"];
	Particle_Ptr const& D1920m   = P["D1920m"];
	Particle_Ptr const& N1440p   = P["N1440p"];
	Particle_Ptr const& N1440n   = P["N1440n"];
	Particle_Ptr const& N1520p   = P["N1520p"];
	Particle_Ptr const& N1520n   = P["N1520n"];
	Particle_Ptr const& Proton   = P["proton"];
	Particle_Ptr const& Neutron  = P["neutron"];
	Particle_Ptr const& Pi_Plus  = P["pi+"];
	Particle_Ptr const& Pi_Minus = P["pi-"];
	Particle_Ptr const& Pi_Zero  = P["pi0"];

	std::vector<Particle_Ptr> resonances = {D1232pp, D1232p, D1232n, D1232m, D1600pp, D1600p, D1600n, D1600m, N1440p, N1440n, N1520p, N1520n, D1920pp, D1920p, D1920n, D1920m};

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
	else if( form_factor == "breit_wigner" ){
		ff = breit_wigner;
	}
	else{
		critical_error("Unknown form factor");
	}

	for( auto& resonance : resonances ){
		resonance->user_data("form_factor", ff);
	}

	if( !cmd.is_enabled("D1232") &&  !cmd.is_enabled("N1440") && !cmd.is_enabled("N1520") && !cmd.is_enabled("N1535") ){
		critical_error("No particle activated.\n");
	}

	init_vertices(P);

	if( cmd.as_string("process") == CMD_PROCESS_ELASTIC_SCATTERING ){
		std::vector<Feynman_Diagram_Ptr> diagrams_proton_pi_plus;
		//std::vector<double> values = {1.321,1.362,1.390,1.417,1.427,1.443,1.450,1.462,1.470,1.481,1.495,1.501,1.512,1.528,1.540,1.561,1.562,1.572,1.586,1.612,1.621,1.638,1.641,1.643,1.669,1.673,1.688,1.694,1.716,1.738,1.769,1.777,1.783,1.791,1.821,1.851,1.839,1.878,1.881,1.896,1.903,1.911,1.927,1.945,1.955,1.969,1.978,2.016,2.020,2.057,2.071,2.089,2.102,2.115,2.155,2.189,2.194,2.206,2.240,2.286,2.306,2.520,2.556,2.575,2.778,2.785,2.868,2.900,3.207,3.487,3.696,3.999,4.105,4.173,4.601,4.781,4.916,4.992,5.355,5.561,5.678,5.968,9.733,11.500,13.732,16.236,18.147,19.396,21.680};
		std::vector<Particle_Ptr> particles;
		if( cmd.is_enabled("D1232") ){
			particles.push_back(D1232pp);
		}
		if( cmd.is_enabled("D1600") ){
			particles.push_back(D1600pp);
		}
		if( cmd.is_enabled("D1920") ){
			particles.push_back(D1920pp);
		}
		for( auto const& particle : particles ){
			if( s_channel ){
				if( particle->spin().j() == 1.5 ){
					particle->width([&](double p2){
						auto const m = std::sqrt(p2);
						auto const ff = particle->user_data<FORM_FACTOR_FUNCTION>("form_factor")(particle, Proton, Pi_Plus, m);
						return particle->width() * ff * ff * dyson_factor_32p(particle, Proton, Pi_Plus, m);
					});
				}
				diagrams_proton_pi_plus.push_back(
						create_diagram(FORMAT("{} s", particle->name()), Scattering_Horizontal_2_to_2, VMP,
						               {Proton, Pi_Plus},
						               {particle},
						               {Pi_Plus, Proton}
						));
			}
		}

		status(FORMAT("2mubarn: {}", 1._2mubarn));
		status(FORMAT("2mbarn: {}", 1._2mbarn));
		status(FORMAT("2barn: {}", 1._2barn));

		Feynman_Process scattering_proton_pi_plus(diagrams_proton_pi_plus);
		scattering_proton_pi_plus.conversion_factor(1._2mbarn);

		status("Scattering proton pi_plus -> proton pi_plus");
		double start = cmd.is_enabled("start") ? cmd.as_double("start") : 1.1;
		double end = cmd.is_enabled("end") ? cmd.as_double("end") : 2.0;
		std::size_t steps = cmd.is_enabled("steps") ? static_cast<std::size_t>(cmd.as_int("steps")) : 100ULL;
		scattering_proton_pi_plus.print_sigma_table(std::cout, start, end, steps);
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
						create_diagram("N1440 s", Scattering_Horizontal_2_to_2, VMP,
						               {Proton, QED::Photon},
						               {N1440p},
						               {Pi_Zero, Proton}
						));
				diagrams_proton_pi_minus.push_back(
						create_diagram("N1440 s", Scattering_Horizontal_2_to_2, VMP,
						               {Neutron, QED::Photon},
						               {N1440n},
						               {Pi_Minus, Proton}
						));
				diagrams_neutron_pi_plus.push_back(
						create_diagram("N1440 s", Scattering_Horizontal_2_to_2, VMP,
						               {Proton, QED::Photon},
						               {N1440p},
						               {Pi_Plus, Neutron}
						));
			}
			if( u_channel ){
				diagrams_proton_pi_zero.push_back(
						create_diagram("N1440 u", Scattering_Vertical_2_to_2, VMP,
						               {Proton, QED::Photon},
						               {N1440p},
						               {Pi_Zero, Proton}
						));
				diagrams_proton_pi_minus.push_back(
						create_diagram("N1440 u", Scattering_Vertical_2_to_2, VMP,
						               {Neutron, QED::Photon},
						               {N1440p},
						               {Pi_Minus, Proton}
						));
				diagrams_neutron_pi_plus.push_back(
						create_diagram("N1440 u", Scattering_Vertical_2_to_2, VMP,
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
						create_diagram("N1520 s", Scattering_Horizontal_2_to_2, VMP,
						               {Proton, QED::Photon},
						               {N1520p},
						               {Pi_Zero, Proton}
						));
				diagrams_proton_pi_minus.push_back(
						create_diagram("N1520 s", Scattering_Horizontal_2_to_2, VMP,
						               {Neutron, QED::Photon},
						               {N1520n},
						               {Pi_Minus, Proton}
						));
				diagrams_neutron_pi_plus.push_back(
						create_diagram("N1520 s", Scattering_Horizontal_2_to_2, VMP,
						               {Proton, QED::Photon},
						               {N1520p},
						               {Pi_Plus, Neutron}
						));
			}
			if( u_channel ){
				diagrams_proton_pi_zero.push_back(
						create_diagram("N1520 u", Scattering_Vertical_2_to_2, VMP,
						               {Proton, QED::Photon},
						               {N1520p},
						               {Pi_Zero, Proton}
						));
				diagrams_proton_pi_minus.push_back(
						create_diagram("N1520 u", Scattering_Vertical_2_to_2, VMP,
						               {Neutron, QED::Photon},
						               {N1520p},
						               {Pi_Minus, Proton}
						));
				diagrams_neutron_pi_plus.push_back(
						create_diagram("N1520 u", Scattering_Vertical_2_to_2, VMP,
						               {Proton, QED::Photon},
						               {N1520n},
						               {Pi_Plus, Neutron}
						));
			}
		}

		Feynman_Process scattering_proton_pi_zero(diagrams_proton_pi_zero);
		Feynman_Process scattering_proton_pi_minus(diagrams_proton_pi_minus);
		Feynman_Process scattering_neutron_pi_plus(diagrams_neutron_pi_plus);

		Timer stopwatch;
		stopwatch.start();

		double sqrt_s = cmd.as_double("sqrt_s");

		std::vector<double> const values = {1.175, 1.2,1.225,1.25,1.275,1.3,1.325,1.35,1.375,1.4,1.41,1.42,1.43,1.44,1.45,1.46,1.47,1.48,1.49,1.5,1.51,1.52,1.53,1.54,1.55,1.56,1.57,1.58,1.59,1.6,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2};
		status("Scattering proton gamma -> proton pi0");
		status(FORMAT("sqrt_s = {}", sqrt_s));
//		scattering_proton_pi_zero.print_dsigma_dcos_table(std::cout, sqrt_s, 0.1);
		scattering_proton_pi_zero.print_sigma_table(std::cout, values);
		std::cout << "\n\n";
		status("Scattering neutron gamma -> proton pi-");
		status(FORMAT("sqrt_s = {}", sqrt_s));
//		scattering_proton_pi_minus.print_dsigma_dcos_table(std::cout, sqrt_s, 0.1);
		scattering_proton_pi_minus.print_sigma_table(std::cout, values);
		status("Scattering proton gamma -> neutron pi+");
		status(FORMAT("sqrt_s = {}", sqrt_s));
//		scattering_neutron_pi_plus.print_dsigma_dcos_table(std::cout, sqrt_s, 0.1);
		scattering_neutron_pi_plus.print_sigma_table(std::cout, {{1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0}});
		scattering_neutron_pi_plus.print_sigma_table(std::cout, values);
		stopwatch.stop();
		std::cout << "Time: " << stopwatch.time<std::chrono::seconds>() << "\n";
	}
	return EXIT_SUCCESS;
}