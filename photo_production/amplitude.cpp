#include "effective_lagrangian_model.hpp"
#include "form_factors.hpp"

#include <feynumeric/feynumeric.hpp>

#include <iostream>

int main(int argc, char** argv){
	using namespace Feynumeric;

	Command_Line_Manager cmd(argc, argv);
	cmd.register_command("particle_file", true, "file with particle parameters");
	cmd.register_command("coupling_constants", true, "file with coupling constants");
	cmd.register_command("start", std::string("1.0"), "starting point");
	cmd.register_command("end", std::string("3.0"), "end value");
	cmd.register_command("steps", std::string("200"), "steps");
	cmd.register_command("form_factor", Form_Factor::CMD_FORM_FACTOR_NONE, FORMAT("which form factor to use ({}, {}, {}, {}, {}, {})", Form_Factor::CMD_FORM_FACTOR_NONE, Form_Factor::CMD_FORM_FACTOR_CASSING, Form_Factor::CMD_FORM_FACTOR_CUTKOSKY, Form_Factor::CMD_FORM_FACTOR_MANLEY, Form_Factor::CMD_FORM_FACTOR_MONIZ, Form_Factor::CMD_FORM_FACTOR_BREIT_WIGNER));
	cmd.register_command("particle", true);
	cmd.crash_on_missing_mandatory_command();

	Particle_Manager P(cmd.as_string("particle_file"));
	auto const& Proton = P.get("proton");
	auto const& Neutron = P.get("neutron");
	auto const& Pi_Zero = P.get("pi0");
	auto const& Pi_Minus = P.get("pi-");
	auto const& Pi_Plus = P.get("pi+");
	init_vertices(P, cmd.as_string("coupling_constants"));


	std::string particle_name = cmd.as_string("particle");

	Particle_Ptr particle = nullptr;

	Feynman_Diagram_Ptr scattering_diagram = nullptr;
	Feynman_Diagram_Ptr decay_diagram = nullptr;

	std::vector<double> s_values;

	if( particle_name == "Rho" ){
		s_values = lin_space(1.1, 2.1, 17);
		particle = P.get("rho0");
//        auto const coupl_str = coupling_string(particle->name(), "N", "N");
//        couplings.set(coupl_str, 1.);
		scattering_diagram = create_diagram(FORMAT("{} u", particle->name()), u_channel, VMP,
		                                    {Proton, Pi_Plus},
		                                    {particle},
		                                    {Pi_Plus, Proton}
		);
	}
    else if( particle_name == "f0" ){
            s_values = lin_space(1.1, 2.1, 17);
            particle = P.get("f0_500");
            auto const coupl_str = coupling_string(particle->name(), "N", "N");
            couplings.set(coupl_str, 1.);
            scattering_diagram = create_diagram(FORMAT("{} u", particle->name()), u_channel, VMP,
                                                {Proton, Pi_Plus},
                                                {particle},
                                                {Pi_Plus, Proton}
            );
        }
    else if( particle_name == "photon" ){
        s_values = lin_space(1.08, 2.08, 100);
        particle = QED::Photon;
        auto const coupl_str = coupling_string(particle->name(), "N", "N");
        couplings.set(coupl_str, 1.);
        scattering_diagram = create_diagram(FORMAT("{} u", particle->name()), u_channel, VMP,
                                            {Proton, Pi_Plus},
                                            {particle},
                                            {Pi_Plus, Proton}
        );
    }
	else{
		particle = P.get(FORMAT("{}pp", cmd.as_string("particle")));
		s_values = weighted_space(1.1, particle->mass() - particle->width(), particle->mass() + particle->width(), 2.1, 17);
        auto const coupl_str = coupling_string(particle->name(), "N", "Pion");
        couplings.set(coupl_str, 1.);
		scattering_diagram = create_diagram(FORMAT("{} s", particle->name()), s_channel, VMP,
		                                    {Proton, Pi_Plus},
		                                    {particle},
		                                    {Pi_Plus, Proton}
		);
		decay_diagram = create_diagram(FORMAT("{} decay", particle->name()), Decay_1_to_2, VMP,
		                               {particle},
		                               {},
		                               {Pi_Plus, Proton});
	}


	particle->user_data("form_factor", Form_Factor::identity);

	if( decay_diagram ){
		Feynman_Process process_decay({decay_diagram});
		auto result_decay = process_decay.decay_amplitude(s_values, 4);
		for( std::size_t i = 0; i < result_decay.size(); ++i ){
			auto const& row = result_decay[i];
			for( std::size_t j = 0; j < row.size(); ++j ){
				auto const& poly = row[j];
				poly.save(FORMAT("data/polynomials/polynomial_decay_{}_{}_{}.txt", particle->name(), i, j));
			}
		}
	}

	if( scattering_diagram ){
		Feynman_Process process_scattering({scattering_diagram});

        std::vector<std::vector<Polynomial>> result_scattering;

        if( particle_name == "photon"){
            result_scattering = process_scattering.scattering_amplitude(s_values, {40, 8});
        }
        else{
            result_scattering = process_scattering.scattering_amplitude(s_values, {6, 8});
        }

		for( std::size_t i = 0; i < result_scattering.size(); ++i ){
			auto const& row = result_scattering[i];
			for( std::size_t j = 0; j < row.size(); ++j ){
				auto const& poly = row[j];
				poly.save(FORMAT("data/polynomials/polynomial_scattering_{}_{}_{}.txt", particle->name(), i, j));
                if( i == 0 ){
                    std::cout << poly.to_string('s') << "\n";
                }
                if( i == 1 ){
                    std::cout << poly.to_string('c') << "\n";
                }

			}
		}
	}


}