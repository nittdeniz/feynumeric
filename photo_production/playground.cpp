#include <feynumeric/feynumeric.hpp>

#include "couplings.hpp"
#include "effective_lagrangian_model.hpp"
#include "form_factors.hpp"

#include <fstream>
#include <sstream>
#include <string>

int main(int argc, char** argv){
	using namespace Feynumeric;
	using namespace Feynumeric::Units;
	using namespace std::placeholders;

	Command_Line_Manager cmd(argc, argv);
	cmd.register_command("particle_file", true, "file with particle parameters");
	cmd.register_command("coupling_constants", true, "file with coupling constants");
	cmd.register_command("start", std::string("1.0"), "starting point");
	cmd.register_command("end", std::string("3.0"), "end value");
	cmd.register_command("steps", std::string("200"), "steps");
	cmd.register_command("particle", true);

	cmd.crash_on_missing_mandatory_command();

	Particle_Manager P(cmd.as_string("particle_file"));
	auto const& Proton = P.get("proton");
	auto const& Neutron = P.get("neutron");
	auto const& Pi_Zero = P.get("pi0");
	auto const& Pi_Minus = P.get("pi-");
	auto const& Pi_Plus = P.get("pi+");

	auto particle = P.get(FORMAT("{}pp", cmd.as_string("particle")));

	init_vertices(P, cmd.as_string("coupling_constants"));

	double const start = cmd.as_double("start");
	double const end   = cmd.as_double("end");
	int    const steps = cmd.as_int("steps");

//	particle->feynman_virtual = [](std::shared_ptr<Graph_Edge> e, Kinematics const& kin){ return Projector(e, kin); };
	particle->user_data("form_factor", identity);

	init_vertices(P, cmd.as_string("coupling_constants"));
	Timer stopwatch;
	stopwatch.start();

	std::vector<double> cos_theta_values;
	cos_theta_values.reserve(201);
	cos_theta_values.push_back(-0.9999);
	for( std::size_t i = 1; i < 200; ++i ){
		cos_theta_values.push_back(-1 + 0.01 * i);
	}
	cos_theta_values.push_back(0.9999);

	std::stringstream output_proc1;
	std::stringstream output_proc2;
	std::stringstream output_proc3;
	std::stringstream output_proc4;
	std::stringstream output_proc5;

	auto coupl_str = coupling_string(cmd.as_string("particle"), "N", "Pion");
	couplings.set(coupl_str, 1.);

	for( int i = 0; i < 1; ++i ){
		double sqrt_s = 1.23;
		if( cmd.as_string("particle") == "D1232" ){
			auto dummy = std::make_shared<Particle>(*particle);
			auto scattering = create_diagram(FORMAT("{} s", particle->name()), s_channel, VMP,
			                                 {Proton, Pi_Plus},
			                                 {particle},
			                                 {Pi_Plus, Proton}
			);
			Feynman_Process process({scattering});
			auto result = process.dsigma_dcos_table(sqrt_s, 200ULL);
			for( auto [key, array] : result ){
				std::cout << key;
				std::cout << "\t" << array[0];
				std::cout << "\n";
			}
			std::cout << "integrated: \n";
			process.print_sigma_table(std::cout, {sqrt_s}, 0.01);
		}
	}
}

