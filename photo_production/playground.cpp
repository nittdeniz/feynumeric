#include "effective_lagrangian_model.hpp"
#include "form_factors.hpp"


#include <feynumeric/feynumeric.hpp>
#include <feynumeric/utility.hpp>

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
	cmd.register_command("particle", true);
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

	auto particle = P.get(FORMAT("{}pp", cmd.as_string("particle")));

	Polynomial test({3.+5.i});
	func_t<1> ffff = test;
	FPolynomial<1> ftest(std::vector<func_t<1>>{{ffff}});
	std::function<Complex(Complex)> cj = [](Complex c){ return std::conj(c);};
	auto ftestc = ftest >> cj;

//	std::cout << "@@1: " << ftest(1., 1.) << "\n";
//	std::cout << "@@2: " << ftestc(1., 1.) << "\n";

	/* width */

	auto n_spin_states = particle->spin().n_states() * Proton->spin().n_states();
	std::vector<std::vector<Polynomial>> width_polynomials(1);
	for( std::size_t i = 0; i < n_spin_states; ++i ){
		Polynomial temp;
		temp.load(FORMAT("data/polynomials/polynomial_decay_{}_0_{}.txt", particle->name(), i));
		width_polynomials[0].push_back(temp);
	}

	Amplitude<0> M_width(width_polynomials, {particle}, {}, {Proton, Pi_Plus});

//	std::vector<double> wv;
	auto sqrt_s_values = lin_space(1.1, 2.1, 100);
//	auto width_result = M_width.width(sqrt_s_values);
//	for( auto const& [s, w] : width_result ){
//		std::cout << s << " " << w << "\n";
//	}


	/* scattering */
	n_spin_states = Proton->spin().n_states() * Proton->spin().n_states();
	std::vector<std::vector<Polynomial>> scattering_polynomials(2);
	for( std::size_t i = 0; i < n_spin_states; ++i ){
		Polynomial temp_s, temp_c;
		temp_s.load(FORMAT("data/polynomials/polynomial_scattering_{}_0_{}.txt", particle->name(), i));
		temp_c.load(FORMAT("data/polynomials/polynomial_scattering_{}_1_{}.txt", particle->name(), i));
		scattering_polynomials[0].push_back(temp_s);
//		std::cout << temp_s.to_string('s') << "\n";
		scattering_polynomials[1].push_back(temp_c);
//		std::cout << temp_c.to_string('c') << "\n";
	}

//	std::cout << scattering_polynomials[0][1].to_string('s') << "\n";
//	std::cout << scattering_polynomials[1][1].to_string('c') << "\n";

	func_t<1> phase_space = [&](double sqrt_s){
		auto qout = momentum(sqrt_s, Proton->mass(), Pi_Plus->mass());
		auto qin = qout;
		return phase_space2(4, sqrt_s, qout, qin);
	};

	func_t<1> breit_wigner = [&](double sqrt_s){
		return 1./(sqrt_s*sqrt_s - particle->mass() * particle->mass() + 1.i * sqrt_s * M_width.width({sqrt_s})[sqrt_s]);
//		return 1./(sqrt_s*sqrt_s - particle->mass() * particle->mass());
	};

	Amplitude<1> M_scattering(scattering_polynomials, {Proton, Pi_Plus}, {particle},{Proton, Pi_Plus});
	//M_scattering.scale(phase_space);

	func_t<1> form_factor = [&](double sqrt_s){
		auto const lambda = 0.8;
		auto const l4 = std::pow(lambda, 4);
		auto const delta = sqrt_s * sqrt_s - particle->mass() * particle->mass();
		return l4/(delta*delta + l4);
	};

	auto add = M_scattering + M_scattering;

//	M_scattering.scale(breit_wigner);
//	M_scattering.scale(form_factor);

//	add.scale(breit_wigner);
//	add.scale(form_factor);

	double sqrt_s = 1.2;
	double cos    = 0;


	auto result_scattering = M_scattering.scattering(sqrt_s_values);

	auto result_add = add.scattering(sqrt_s_values);



	using namespace Feynumeric::Units;
	std::cout << "{";
	for( auto const& [x, y] : result_scattering ){
		std::cout << "{" << x << "," << y * 1._2mbarn << "},";
	}
	std::cout << "\b};\n";

	std::cout << "{";
	for( auto const& [x, y] : result_add ){
		std::cout << "{" << x << "," << y * 1._2mbarn << "},";
	}
	std::cout << "\b};";

}