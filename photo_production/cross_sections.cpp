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
	cmd.register_command("form_factor", std::string("none"), "which form factor to use (none, cassing, manley, monitz)");
	cmd.register_command("channel", std::string("s"), "which channel to use [s, t, u, c] or any combination.");
	cmd.register_command("sqrt_s", true, "the energy in the center of mass frame in GeV");

	cmd.crash_on_missing_mandatory_command();

	double sqrt_s = cmd.as_double("sqrt_s");

	auto const& channel = cmd.as_string("channel");
	bool s_channel = channel.find("s") != std::string::npos;
	bool u_channel = channel.find("u") != std::string::npos;
	bool t_channel = channel.find("t") != std::string::npos;
	bool c_channel = channel.find("c") != std::string::npos;




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
					create_diagram("N1440 s_channel1", Double_Wrench, VMP,
					               {Proton, QED::Photon},
					               {N1440p},
					               {Pi_Zero, Proton}
					));
			diagrams_proton_pi_minus.push_back(
					create_diagram("N1440 s_channel2", Double_Wrench, VMP,
					               {Neutron, QED::Photon},
					               {N1440n},
					               {Pi_Minus, Proton}
					));
			diagrams_neutron_pi_plus.push_back(
					create_diagram("N1440 s_channel3", Double_Wrench, VMP,
					               {Proton, QED::Photon},
					               {N1440p},
					               {Pi_Plus, Neutron}
					));
		}
		if( u_channel ){
			diagrams_proton_pi_zero.push_back(
					create_diagram("N1440 u_channel1", X_Man, VMP,
					               {Proton, QED::Photon},
					               {N1440p},
					               {Pi_Zero, Proton}
					));
			diagrams_proton_pi_minus.push_back(
					create_diagram("N1440 u_channel2", X_Man, VMP,
					               {Neutron, QED::Photon},
					               {N1440n},
					               {Pi_Minus, Proton}
					));
			diagrams_neutron_pi_plus.push_back(
					create_diagram("N1440 u_channel3", X_Man, VMP,
					               {Proton, QED::Photon},
					               {N1440p},
					               {Pi_Plus, Neutron}
					));
		}
	}
	if( cmd.is_enabled("N1520") )
	{
		if( s_channel ){
			diagrams_proton_pi_zero.push_back(
					create_diagram("N1520 s_channel1", Double_Wrench, VMP,
					               {Proton, QED::Photon},
					               {N1520p},
					               {Pi_Zero, Proton}
					));
			diagrams_proton_pi_minus.push_back(
					create_diagram("N1520 s_channel2", Double_Wrench, VMP,
					               {Neutron, QED::Photon},
					               {N1520n},
					               {Pi_Minus, Proton}
					));
			diagrams_neutron_pi_plus.push_back(
					create_diagram("N1520 s_channel3", Double_Wrench, VMP,
					               {Proton, QED::Photon},
					               {N1520p},
					               {Pi_Plus, Neutron}
					));
		}
		if( u_channel ){
			diagrams_proton_pi_zero.push_back(
					create_diagram("N1520 u_channel1", X_Man, VMP,
					               {Proton, QED::Photon},
					               {N1520p},
					               {Pi_Zero, Proton}
					));
			diagrams_proton_pi_minus.push_back(
					create_diagram("N1520 u_channel2", X_Man, VMP,
					               {Neutron, QED::Photon},
					               {N1520n},
					               {Pi_Minus, Proton}
					));
			diagrams_neutron_pi_plus.push_back(
					create_diagram("N1520 u_channel3", X_Man, VMP,
					               {Proton, QED::Photon},
					               {N1520p},
					               {Pi_Plus, Neutron}
					));
		}
	}

	Feynman_Process scattering_proton_pi_zero(diagrams_proton_pi_zero);
	Feynman_Process scattering_proton_pi_minus(diagrams_proton_pi_minus);
	Feynman_Process scattering_neutron_pi_plus(diagrams_neutron_pi_plus);

	std::cout << "Scattering proton gamma -> proton pi0\n";
	std::cout << "sqrt_s = " << sqrt_s << "\n";
	scattering_proton_pi_zero.print_dsigma_dcos_table(std::cout, sqrt_s, 0.1);
	std::cout << "Scattering neutron gamma -> proton pi-\n";
	std::cout << "sqrt_s = " << sqrt_s << "\n";
	scattering_proton_pi_minus.print_dsigma_dcos_table(std::cout, sqrt_s, 0.1);
	std::cout << "Scattering proton gamma -> neutron pi+\n";
	std::cout << "sqrt_s = " << sqrt_s << "\n";
	scattering_proton_pi_zero.print_dsigma_dcos_table(std::cout, sqrt_s, 0.1);
	return EXIT_SUCCESS;
}