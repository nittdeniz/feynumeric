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

	Command_Line_Manager cmd(argc, argv);

	cmd.register_command("particle_file", true, "file with particle parameters");
	cmd.register_command("form_factor", std::string("none"), "which form factor to use (none, cassing, manley, moniz)");
	cmd.register_command("channel", std::string("s"), "which channel to use [s, t, u, c] or any combination.");
	cmd.register_command("sqrt_s", true, "the energy in the center of mass frame in GeV");

	cmd.crash_on_missing_mandatory_command();

	double sqrt_s = cmd.as_double("sqrt_s");

	auto const& channel = cmd.as_string("channel");
	bool const s_channel = channel.find('s') != std::string::npos;
	bool const u_channel = channel.find('u') != std::string::npos;
	bool const t_channel = channel.find('t') != std::string::npos;
	bool const c_channel = channel.find('c') != std::string::npos;




	Particle_Manager P(cmd.as_string("particle_file"));
	Particle_Ptr const& N1440p   = P.get("N1440p");
	Particle_Ptr const& N1440n   = P.get("N1440n");
	Particle_Ptr const& N1520p   = P.get("N1520p");
	Particle_Ptr const& N1520n   = P.get("N1520n");
	Particle_Ptr const& Proton   = P.get("proton");
	Particle_Ptr const& Neutron  = P.get("neutron");
	Particle_Ptr const& Pi_Plus  = P.get("pi+");
	Particle_Ptr const& Pi_Minus = P.get("pi-");
	Particle_Ptr const& Pi_Zero  = P.get("pi0");

	std::string const& form_factor = cmd.as_string("form_factor");

	if( form_factor == "none" ){
		N1440p->user_data("form_factor", identity);
		N1440n->user_data("form_factor", identity);
		N1520p->user_data("form_factor", identity);
		N1520n->user_data("form_factor", identity);
	}
	else if( form_factor == "moniz" ){
		N1440p->user_data("form_factor", moniz);
		N1440n->user_data("form_factor", moniz);
		N1520p->user_data("form_factor", moniz);
		N1520n->user_data("form_factor", moniz);
	}
	else if( form_factor == "manley" ){
		N1440p->user_data("form_factor", manley);
		N1440n->user_data("form_factor", manley);
		N1520p->user_data("form_factor", manley);
		N1520n->user_data("form_factor", manley);
	}
	else if( form_factor == "cutkosky" ){
		N1440p->user_data("form_factor", cutkosky);
		N1440n->user_data("form_factor", cutkosky);
		N1520p->user_data("form_factor", cutkosky);
		N1520n->user_data("form_factor", cutkosky);
	}
	else if( form_factor == "cassing" ){
		N1440p->user_data("form_factor", cassing);
		N1440n->user_data("form_factor", cassing);
		N1520p->user_data("form_factor", cassing);
		N1520n->user_data("form_factor", cassing);
	}
	else{
		critical_error("Unknown form factor");
	}

	init_vertices(P);

	// Photo Production
	std::vector<Feynman_Diagram_Ptr> diagrams_proton_pi_zero;
	std::vector<Feynman_Diagram_Ptr> diagrams_proton_pi_minus;
	std::vector<Feynman_Diagram_Ptr> diagrams_neutron_pi_plus;

	if( !cmd.is_enabled("N1440") && !cmd.is_enabled("N1520") ){
		critical_error("No particle activated.\n");
	}

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
	std::cout << "\n\n";
	status("Scattering neutron gamma -> proton pi-");
	status(FORMAT("sqrt_s = {}", sqrt_s));
	scattering_proton_pi_minus.print_dsigma_dcos_table(std::cout, sqrt_s, 0.1);
	status("Scattering proton gamma -> neutron pi+");
	status(FORMAT("sqrt_s = {}", sqrt_s));
	scattering_proton_pi_zero.print_dsigma_dcos_table(std::cout, sqrt_s, 0.1);
	return EXIT_SUCCESS;
}