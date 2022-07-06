#include "effective_lagrangian_model.hpp"
#include "form_factors.hpp"


#include <feynumeric/feynumeric.hpp>
#include <feynumeric/utility.hpp>

#include <chrono>
#include <iostream>
#include <feynumeric/phase_space.hpp>

int main(int argc, char** argv){
	using namespace Feynumeric;

	Command_Line_Manager cmd(argc, argv);
	cmd.register_command("particle_file", true, "file with particle parameters");
	cmd.register_command("coupling_constants", true, "file with coupling constants");
	cmd.register_command("start", std::string("1.0"), "starting point");
	cmd.register_command("end", std::string("3.0"), "end value");
	cmd.register_command("steps", std::string("200"), "steps");
	cmd.register_command("form_factor", Form_Factor::CMD_FORM_FACTOR_NONE, FORMAT("which form factor to use ({}, {}, {}, {}, {}, {})", Form_Factor::CMD_FORM_FACTOR_NONE, Form_Factor::CMD_FORM_FACTOR_CASSING, Form_Factor::CMD_FORM_FACTOR_CUTKOSKY, Form_Factor::CMD_FORM_FACTOR_MANLEY, Form_Factor::CMD_FORM_FACTOR_MONIZ, Form_Factor::CMD_FORM_FACTOR_BREIT_WIGNER));
	cmd.crash_on_missing_mandatory_command();

	Particle_Manager P(cmd.as_string("particle_file"));
	auto const& Proton = P.get("proton");
	auto const& Neutron = P.get("neutron");
	auto const& Pi_Zero = P.get("pi0");
	auto const& Pi_Minus = P.get("pi-");
	auto const& Pi_Plus = P.get("pi+");
	init_vertices(P, cmd.as_string("coupling_constants"));

	std::vector<Polynomial> s_poly;
	std::vector<Polynomial> c_poly;

	std::vector<std::string> resonances = {"D1232", "D1600", "D1620", "D1700", "D1750", "D1900", "D1905", "D1910", "D1920", "D1930", "D1940", "D1950"};

	std::vector<Amplitude<0>> width_amplitudes;
	std::vector<Amplitude<1>> scattering_amplitudes;

	Timer stopwatch;
	stopwatch.start();

	for( auto const& resonance : resonances ){
		if( !cmd.is_enabled(resonance)){
			continue;
		}
		auto particle = P.get(FORMAT("{}pp", resonance));
		/* width */

		auto n_spin_states = particle->spin().n_states() * Proton->spin().n_states();
		std::vector<std::vector<Polynomial>> width_polynomials(1);
		for( std::size_t i = 0; i < n_spin_states; ++i ){
			Polynomial temp;
			temp.load(FORMAT("data/polynomials/polynomial_decay_{}_0_{}.txt", particle->name(), i));
			width_polynomials[0].push_back(temp);
		}

		Amplitude<0> M_width(width_polynomials, {particle}, {}, {Proton, Pi_Plus});
		width_amplitudes.push_back(M_width);

		/* scattering */
		n_spin_states = Proton->spin().n_states() * Proton->spin().n_states();
		std::vector<std::vector<Polynomial>> scattering_polynomials(2);
		for( std::size_t i = 0; i < n_spin_states; ++i ){
			Polynomial temp_s, temp_c;
			temp_s.load(FORMAT("data/polynomials/polynomial_scattering_{}_0_{}.txt", particle->name(), i));
			temp_c.load(FORMAT("data/polynomials/polynomial_scattering_{}_1_{}.txt", particle->name(), i));
			scattering_polynomials[0].push_back(temp_s);
			scattering_polynomials[1].push_back(temp_c);
		}

		func_t<1> phase_space = [&](double sqrt_s){
			auto qout = momentum(sqrt_s, Proton->mass(), Pi_Plus->mass());
			auto qin = qout;
			return phase_space2(4, sqrt_s, qout, qin);
		};

		func_t<1> breit_wigner = [=](double sqrt_s) mutable {
			return 1. / ( sqrt_s * sqrt_s - particle->mass() * particle->mass() +
			              1.i * sqrt_s * M_width.width({sqrt_s})[sqrt_s] );
		};

		Amplitude<1> M_scattering(scattering_polynomials, {Proton, Pi_Plus}, {particle}, {Proton, Pi_Plus});

		func_t<1> form_factor = [=](double sqrt_s) mutable{
			auto const lambda = 0.8;
			auto const l4 = std::pow(lambda, 4);
			auto const delta = sqrt_s * sqrt_s - particle->mass() * particle->mass();
			return l4 / ( delta * delta + l4 );
		};

		M_scattering.scale(breit_wigner);
		M_scattering.scale(form_factor);
		scattering_amplitudes.push_back(M_scattering);

	}

	if( scattering_amplitudes.empty() ){
		error("No particle selected.");
		return EXIT_SUCCESS;
	}

	Amplitude<1> interference = scattering_amplitudes[0];
	for( std::size_t i = 1; i < scattering_amplitudes.size(); ++i ){
		interference = interference + scattering_amplitudes[i];
	}

	stopwatch.stop();
	std::cout << "config done: " << stopwatch.time<std::chrono::milliseconds>()/1000. << "\n";
	stopwatch.start();

	auto sqrt_s_values = lin_space(1.1, 2.1, 10);
	auto result = interference.scattering(sqrt_s_values);

	std::cout << "{";
	for( auto& [s,val] : result ){
		std::cout << "{" << s << "," << val << "},";
	}
	std::cout << "\b};\n";

	stopwatch.stop();
	std::cout << "config done: " << stopwatch.time<std::chrono::milliseconds>()/1000. << "\n";




	//auto result_scattering = M_scattering.scattering(sqrt_s_values);
}