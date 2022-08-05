#include <feynumeric/command_line_manager.hpp>
#include <feynumeric/core.hpp>
#include <feynumeric/table.hpp>
#include <feynumeric/qed.hpp>
#include <feynumeric/messages.hpp>
#include <feynumeric/timer.hpp>

#include <filesystem>

#include "effective_lagrangian_model.hpp"
#include "form_factors.hpp"
#include <omp.h>


#include <iostream>
#include <iomanip>

std::string const CMD_PROCESS_ELASTIC_SCATTERING = "elasticscattering";
std::string const CMD_PROCESS_PHOTO_PRODUCTION   = "photoproduction";
std::string const CMD_CHANNEL_S                  = "s";
std::string const CMD_CHANNEL_T                  = "t";
std::string const CMD_CHANNEL_U                  = "u";
std::string const CMD_CHANNEL_C                  = "c";
std::string const CMD_CROSS_SECTION_TOTAL        = "total";
std::string const CMD_CROSS_SECTION_DIFFERENTIAL = "differential";

int main(int argc, char** argv)
{
	using namespace Feynumeric;
	using namespace Feynumeric::Units;

	Command_Line_Manager cmd(argc, argv);

	cmd.register_command("particle_file", true, "file with particle parameters");
    cmd.register_command("coupling_constants", true, "file with coupling constants");
	cmd.register_command("form_factor", Form_Factor::CMD_FORM_FACTOR_NONE, FORMAT("which form factor to use ({}, {}, {}, {}, {}, {})", Form_Factor::CMD_FORM_FACTOR_NONE, Form_Factor::CMD_FORM_FACTOR_CASSING, Form_Factor::CMD_FORM_FACTOR_CUTKOSKY, Form_Factor::CMD_FORM_FACTOR_MANLEY, Form_Factor::CMD_FORM_FACTOR_MONIZ, Form_Factor::CMD_FORM_FACTOR_BREIT_WIGNER));
	cmd.register_command("channel", CMD_CHANNEL_S, "which channel to use [s, t, u, c] or any combination.");
	cmd.register_command("process", CMD_PROCESS_PHOTO_PRODUCTION, FORMAT("which process to use: {} {}", CMD_PROCESS_PHOTO_PRODUCTION, CMD_PROCESS_ELASTIC_SCATTERING));
	cmd.register_command("sqrt_s", false, "the energy in the center of mass frame in GeV");
	cmd.register_command("help", false, "list all command line parameters");
	cmd.register_command("cross_section", CMD_CROSS_SECTION_TOTAL, FORMAT("get {} or {} cross section", CMD_CROSS_SECTION_DIFFERENTIAL, CMD_CROSS_SECTION_TOTAL));

	if( cmd.is_enabled("help") ){
		return EXIT_SUCCESS;
	}

	cmd.crash_on_missing_mandatory_command();

	auto const& channel = cmd.as_string("channel");
	bool const s_channel_enabled = channel.find('s') != std::string::npos;
	bool const t_channel_enabled = channel.find('t') != std::string::npos;
	bool const u_channel_enabled = channel.find('u') != std::string::npos;
	bool const c_channel_enabled = channel.find('c') != std::string::npos;

	Particle_Manager P(cmd.as_string("particle_file"));
	Particle_Ptr const& Proton   = P.get("proton");
	Particle_Ptr const& Neutron  = P.get("neutron");
	Particle_Ptr const& Pi_Plus  = P.get("pi+");
	Particle_Ptr const& Pi_Minus = P.get("pi-");
	Particle_Ptr const& Pi_Zero  = P.get("pi0");


	std::string const& form_factor = cmd.as_string("form_factor");

	FORM_FACTOR_FUNCTION ff;
	if( Form_Factor::ff_dict.contains(form_factor) ){
	    ff = Form_Factor::ff_dict[form_factor];
	}
	else{
		critical_error("Unknown form factor");
	}

	status(FORMAT("Form factor: {}", form_factor));

	std::vector<std::string> const nucleon_resonances = {"N1440", "N1520", "N1535", "N1650", "N1675", "N1680", "N1700", "N1710", "N1720", "N1875", "N1880", "N1895", "N1900", "N2060", "N2100", "N2120", "N2190", "N2220", "N2250", "N2600"};
	std::vector<std::string> const delta_resonances   = {"D1232", "D1600", "D1620", "D1700", "D1900", "D1905", "D1910", "D1920", "D1930", "D1940", "D1950", "D2200", "D2420"};

	std::vector<std::string> resonances;
	resonances.insert(resonances.end(), nucleon_resonances.cbegin(), nucleon_resonances.cend());
	resonances.insert(resonances.end(), delta_resonances.cbegin(), delta_resonances.cend());

	auto pp_string = [](std::string const& p){ return FORMAT("{}pp", p);};
	auto p_string  = [](std::string const& p){ return FORMAT("{}p", p);};
	auto n_string  = [](std::string const& p){ return FORMAT("{}n", p);};
	auto m_string  = [](std::string const& p){ return FORMAT("{}m", p);};


	for( auto const& nucleon_resonance : nucleon_resonances ){
        std::ifstream ifs(FORMAT("./data/dyson_factors/{}_breit_wigner.txt", nucleon_resonance));
        std::cout << FORMAT("./data/dyson_factors/{}_breit_wigner.txt", nucleon_resonance) << "\n";
		auto const& Np = p_string(nucleon_resonance);
		auto const& Nn = n_string(nucleon_resonance);
		if( P.exists(Np) ){
			P[Np]->user_data("form_factor", ff);
			if( ifs ){
                Table dyson({});
                ifs >> dyson;
                P[Np]->user_data("dyson_factors", dyson);
			}
		}
		if( P.exists(Nn) ){
			P[Nn]->user_data("form_factor", ff);
            if( ifs ){
                Table dyson({});
                ifs >> dyson;
                P[Np]->user_data("dyson_factors", dyson);
            }
		}
	}
	for( auto const& delta_resonance : delta_resonances ){
//	    if( delta_resonance == "D1232" ){
//            ff_temp_str = "breit_wigner";
//            ff = breit_wigner;
//	    }
        std::ifstream ifs(FORMAT("./data/dyson_factors/{}_{}.txt", delta_resonance.substr(0, 5), form_factor));
        std::cout << FORMAT("./data/dyson_factors/{}_{}.txt", delta_resonance.substr(0, 5), form_factor) << "\n";
		auto const& Dpp = pp_string(delta_resonance);
		auto const& Dp  = p_string(delta_resonance);
		auto const& Dn  = n_string(delta_resonance);
		auto const& Dm  = m_string(delta_resonance);
		if( P.exists(Dpp) ){
			P[Dpp]->user_data("form_factor", ff);
            if( ifs ){
                std::cout << "yes\n";
                Table dyson({});
                ifs >> dyson;
                P[Dpp]->user_data("dyson_factors", dyson);
            }
		}
		if( P.exists(Dp) ){
			P[Dp]->user_data("form_factor", ff);
            if( ifs ){
                Table dyson({});
                ifs >> dyson;
                P[Dp]->user_data("dyson_factors", dyson);
            }
		}
		if( P.exists(Dn) ){
			P[Dn]->user_data("form_factor", ff);
            if( ifs ){
                Table dyson({});
                ifs >> dyson;
                P[Dn]->user_data("dyson_factors", dyson);
            }
		}
		if( P.exists(Dm) ){
			P[Dm]->user_data("form_factor", ff);
            if( ifs ){
                Table dyson({});
                ifs >> dyson;
                P[Dm]->user_data("dyson_factors", dyson);
            }
		}
	}

	init_vertices(P, cmd.as_string("coupling_constants"));
	Timer stopwatch;
	stopwatch.start();
	if( cmd.as_string("process") == CMD_PROCESS_ELASTIC_SCATTERING ){
		std::vector<Feynman_Diagram_Ptr> diagrams_proton_pi_plus;
		//std::vector<double> values = {1.321,1.362,1.390,1.417,1.427,1.443,1.450,1.462,1.470,1.481,1.495,1.501,1.512,1.528,1.540,1.561,1.562,1.572,1.586,1.612,1.621,1.638,1.641,1.643,1.669,1.673,1.688,1.694,1.716,1.738,1.769,1.777,1.783,1.791,1.821,1.851,1.839,1.878,1.881,1.896,1.903,1.911,1.927,1.945,1.955,1.969,1.978,2.016,2.020,2.057,2.071,2.089,2.102,2.115,2.155,2.189,2.194,2.206,2.240,2.286,2.306,2.520,2.556,2.575,2.778,2.785,2.868,2.900,3.207,3.487,3.696,3.999,4.105,4.173,4.601,4.781,4.916,4.992,5.355,5.561,5.678,5.968,9.733,11.500,13.732,16.236,18.147,19.396,21.680};
		std::vector<Particle_Ptr> particles;
		for( auto const& R : resonances ){
			if( cmd.is_enabled(R) ){
				status(FORMAT("{} activated", R));
				if( R[0] == 'D' ){
					particles.push_back(P.get(pp_string(R)));
				}
				particles.push_back(P.get(n_string(R)));
			}
		}

		for( auto const& particle : particles ){
		    if( particle->exists_user_data("dyson_factors") ){
                particle->width([&](double p2){
                    auto result = particle->user_data<Table>("dyson_factors").interpolate(std::sqrt(p2));
                    return result;
                });
		    }
			else{
				particle->width([&](double p2){
					return particle->width();
				});
			}
			if( s_channel_enabled && particle->charge() == 2){
				diagrams_proton_pi_plus.push_back(
						create_diagram(FORMAT("{} s", particle->name()), s_channel, VMP,
						               {Proton, Pi_Plus},
						               {particle},
						               {Pi_Plus, Proton}
						));
			}
			if( t_channel_enabled && particle->charge() == 0){
                diagrams_proton_pi_plus.push_back(
                        create_diagram(FORMAT("{} t", particle->name()), t_channel, VMP,
                                       {Proton, Pi_Plus},
                                       {particle},
                                       {Pi_Plus, Proton}
                        )
                );
			}
		}


		if( t_channel_enabled && cmd.exists("Nucleon") ){
			diagrams_proton_pi_plus.push_back(
					create_diagram(FORMAT("{} u", Neutron->name()), t_channel, VMP,
					               {Proton, Pi_Plus},
					               {Neutron},
					               {Pi_Plus, Proton}
					)
			);
		}
        if( u_channel_enabled ){
            diagrams_proton_pi_plus.push_back(
                    create_diagram(FORMAT("{} u", "rho0"), u_channel, VMP,
                                   {Proton, Pi_Plus},
                                   {P.get("rho0")},
                                   {Pi_Plus, Proton}
                                   )
            );
            diagrams_proton_pi_plus.push_back(
                    create_diagram(FORMAT("{} u", "f0_500"), u_channel, VMP,
                                   {Proton, Pi_Plus},
                                   {P.get("f0_500")},
                                   {Pi_Plus, Proton}
                    )
            );
        }

		Feynman_Process scattering_proton_pi_plus(diagrams_proton_pi_plus);
		scattering_proton_pi_plus.conversion_factor(1._2mbarn);

		status("Scattering proton pi_plus -> proton pi_plus");
		double start = cmd.exists("start") ? cmd.as_double("start") : 1.1;
		double end = cmd.exists("end") ? cmd.as_double("end") : 2.0;
		std::cout << std::flush;
		std::size_t steps = cmd.exists("steps") ? static_cast<std::size_t>(cmd.as_int("steps")) : 100ULL;
		scattering_proton_pi_plus.print_sigma_table(std::cout, start, end, steps);
//		std::cout << "\n\n";
//        scattering_proton_pi_plus.print_dsigma_dcos_table(std::cout, start, 100ULL);
//        std::cout << "\n\n";
//        scattering_proton_pi_plus.print_dsigma_dcos_table(std::cout, end, 100ULL);

	}
	else if( cmd.as_string("process") == CMD_PROCESS_PHOTO_PRODUCTION ){
		std::vector<Feynman_Diagram_Ptr> diagrams_proton_pi_zero;
		std::vector<Feynman_Diagram_Ptr> diagrams_proton_pi_minus;
		std::vector<Feynman_Diagram_Ptr> diagrams_neutron_pi_plus;


		std::vector<Particle_Ptr> N_plus = {};
		std::vector<Particle_Ptr> N_null = {};


		for( auto const& particle : resonances ){
			if( cmd.is_enabled(particle) ){
				auto const& Rp = P[p_string(particle)];
				auto const& Rn = P[n_string(particle)];
				if( s_channel_enabled ){
					diagrams_proton_pi_zero.push_back(
							create_diagram(FORMAT("{} s", particle), s_channel, VMP,
							               {Proton, QED::Photon},
							               {Rp},
							               {Pi_Zero, Proton}
							));
					diagrams_proton_pi_minus.push_back(
							create_diagram(FORMAT("{} s", particle), s_channel, VMP,
							               {Neutron, QED::Photon},
							               {Rn},
							               {Pi_Minus, Proton}
							));
					diagrams_neutron_pi_plus.push_back(
							create_diagram(FORMAT("{} s", particle), s_channel, VMP,
							               {Proton, QED::Photon},
							               {Rp},
							               {Pi_Plus, Neutron}
							));
				}
				if( t_channel_enabled ){
					diagrams_proton_pi_zero.push_back(
							create_diagram(FORMAT("{} u", particle), s_channel, VMP,
							               {Proton, QED::Photon},
							               {Rp},
							               {Pi_Zero, Proton}
							));
					diagrams_proton_pi_minus.push_back(
							create_diagram(FORMAT("{} u", particle), s_channel, VMP,
							               {Neutron, QED::Photon},
							               {Rp},
							               {Pi_Minus, Proton}
							));
					diagrams_neutron_pi_plus.push_back(
							create_diagram(FORMAT("{} u", particle), s_channel, VMP,
							               {Proton, QED::Photon},
							               {Rn},
							               {Pi_Plus, Neutron}
							));
				}
			}
		}

		Feynman_Process scattering_proton_pi_zero(diagrams_proton_pi_zero);
		Feynman_Process scattering_proton_pi_minus(diagrams_proton_pi_minus);
		Feynman_Process scattering_neutron_pi_plus(diagrams_neutron_pi_plus);

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
	}
	stopwatch.stop();
	std::cout << "Time: " << std::setw(5) <<  stopwatch.time<std::chrono::milliseconds>()/1000. << "s\n";
	return EXIT_SUCCESS;
}