#include "effective_lagrangian_model.hpp"
#include "form_factors.hpp"

#include <feynumeric/matrix.hpp>
#include <feynumeric/polynomial.hpp>
#include <feynumeric/feynumeric.hpp>
#include <feynumeric/momentum.hpp>

#include <iomanip>
#include <iostream>

int main(int argc, char** argv){
	using namespace Feynumeric;

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

	Timer stopwatch;
	stopwatch.start();


	auto particle = P.get(FORMAT("{}pp", particle_name));
	particle->user_data("form_factor", identity);
	particle->feynman_virtual = [](std::shared_ptr<Graph_Edge> e, Kinematics const& kin){ return Projector(e, kin); };

	auto dummy = std::make_shared<Particle>(*particle);

	auto decay_diagram = create_diagram(FORMAT("{} -> Npi", dummy->name()), Decay_1_to_2, VMP,
	                                    {dummy},{},{Pi_Plus, Proton});

	Feynman_Process decay_process({decay_diagram});

	auto width_amplitudes = decay_process.decay_M_polynomial(dummy, 1.08, 2.08, 4);
	auto conj = width_amplitudes;
	for( auto& M : conj ){
		M = M.conjugate();
	}
	std::vector<Polynomial> width_results;
	for( std::size_t i = 0; i < width_amplitudes.size(); ++i ){
		width_results.push_back(conj[i] * width_amplitudes[i]);
	}
	Polynomial width_result(0);
	for( auto const& wr : width_results ){
		width_result += wr;
	}

	auto values = lin_space(1.08, 2.08, 100);
	values.push_back(1.232);
	std::cout << "{";
	for( auto const& val : values ){
		auto q = momentum(val, Proton->mass(), Pi_Plus->mass());
		auto res = width_result(val).real() * 1./4. * q / ( 8 * M_PI * val * val);
		std::cout << FORMAT("{{{}, {}}},", val, res);
	}
	std::cout << "\b}";

	auto diagram = create_diagram(FORMAT("{} s", particle->name()), s_channel, VMP,
			               {Proton, Pi_Plus},
			               {particle},
			               {Pi_Plus, Proton}
	);

	Feynman_Process process({diagram});

	auto amplitude_list = process.M(1.08,2.1,particle->mass(), 8,6);
	stopwatch.stop();
	std::cout << stopwatch.time<std::chrono::milliseconds>()/1000 << "seconds\n";

	auto conj_amplitude_list = amplitude_list.second;
	for( auto& i : conj_amplitude_list ){
		i = i.conjugate();
	}

	std::vector<Polynomial> multiplied_amplitudes;
	for( std::size_t k = 0; k < amplitude_list.second.size(); ++k ){
		multiplied_amplitudes.
	}

	auto bw = [width_result, particle](double s){
		double const m = particle->mass();
		return 1./(s*s - m*m + 1.i * s * width_result(s));
	};

//	for( auto const& p : list.first ){
//		std::cout << p.to_string('c') << "\n";
//	}
//
//	std::cout << "\n\n@@@@@@@@@@@@@@@@@@@@\n\n";
//
//	for( auto const& p : list.second ){
//		std::cout << p.to_string('s') << "\n";
//	}



}