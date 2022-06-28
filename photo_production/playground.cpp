#include "effective_lagrangian_model.hpp"
#include "form_factors.hpp"

#include <feynumeric/matrix.hpp>
#include <feynumeric/polynomial.hpp>
#include <feynumeric/feynumeric.hpp>

#include <iomanip>
#include <iostream>

int main(int argc, char** argv){
	using namespace Feynumeric;

	Polynomial p(4);
	p.fit({{1.08,-0.610444},{1.205,-5.3766},{1.33,-9.6909},{1.455,-14.6668},{1.58,-20.4498},{1.705,-27.1257},{1.83,-34.7653},{1.955,-43.434},{2.08,-53.1952}});
	std::cout << p.to_string('x') << "\n";

	Command_Line_Manager cmd(argc, argv);
	cmd.register_command("particle_file", true, "file with particle parameters");
	cmd.register_command("coupling_constants", true, "file with coupling constants");
	cmd.register_command("start", std::string("1.0"), "starting point");
	cmd.register_command("end", std::string("3.0"), "end value");
	cmd.register_command("steps", std::string("200"), "steps");
	cmd.register_command("form_factor", CMD_FORM_FACTOR_NONE, FORMAT("which form factor to use ({}, {}, {}, {}, {}, {})", CMD_FORM_FACTOR_NONE, CMD_FORM_FACTOR_CASSING, CMD_FORM_FACTOR_CUTKOSKY, CMD_FORM_FACTOR_MANLEY, CMD_FORM_FACTOR_MONIZ, CMD_FORM_FACTOR_BREIT_WIGNER));

	cmd.crash_on_missing_mandatory_command();

	Particle_Manager P(cmd.as_string("particle_file"));
	auto const& Proton = P.get("proton");
	auto const& Neutron = P.get("neutron");
	auto const& Pi_Zero = P.get("pi0");
	auto const& Pi_Minus = P.get("pi-");
	auto const& Pi_Plus = P.get("pi+");

	init_vertices(P, cmd.as_string("coupling_constants"));

	double const start = cmd.as_double("start");
	double const end   = cmd.as_double("end");
	int    const steps = cmd.as_int("steps");


	std::cout << FORMAT("start: {} end: {} steps: {}\n", start, end, steps);

	std::string particle_name = "D1232";

	auto dummy = std::make_shared<Particle>(*P.get(FORMAT("{}pp", particle_name)));
	dummy->user_data("form_factor", identity);
	auto diagram_pi1 = create_diagram(FORMAT("decay {} to proton pi+", dummy->name()), Decay_1_to_2, VMP,
	                                  {dummy},
	                                  {},
	                                  {Proton, Pi_Plus}
	);
	Feynman_Process process({diagram_pi1});
	auto result = process.decay_M_polynomial(dummy, start, end, 4);

	for( auto& p : result ){
		std::cout << p.to_string('x') << "\n";
	}

}