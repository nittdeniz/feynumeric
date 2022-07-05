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

	auto particle = P.get("D1232pp");

	Polynomial test({3.+5.i});
	func_t<1> ffff = test;
	FPolynomial<1> ftest(std::vector<func_t<1>>{{ffff}});
	std::function<Complex(Complex)> cj = [](Complex c){ return std::conj(c);};
	auto ftestc = ftest >> cj;

	std::cout << "@@1: " << ftest(1., 1.) << "\n";
	std::cout << "@@2: " << ftestc(1., 1.) << "\n";

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
		std::cout << temp_s.to_string('s') << "\n";
		scattering_polynomials[1].push_back(temp_c);
		std::cout << temp_c.to_string('c') << "\n";
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

	M_scattering.scale(breit_wigner);
	//M_scattering.scale(form_factor);

	double sqrt_s = 1.2;
	double cos    = 0;

	auto result_scattering = M_scattering.scattering(sqrt_s_values);



	using namespace Feynumeric::Units;
	std::cout << "{";
	for( auto const& [x, y] : result_scattering ){
		std::cout << "{" << x << "," << y * 1._2mbarn << "},";
	}
	std::cout << "\b};";

//	auto fff = [&](double c){return scattering(c, sqrt_s);};

	//std::cout << "result: " << M_scattering(cos, sqrt_s) << "\n";

	//std::cout << "scattering: " << phase_space(sqrt_s) * integrate(fff, -1, 1) << "\n";


	/*
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


//	Polynomial s({2.0532073895613170, -2.3407191273658303, 1.0085498885716471, -3.2887394432742036, 1.1399531497684108, 0.9129612019750575, 0.3523629278031232});
//	Polynomial c({-19.5810353883248531, 0.0000000000000153, 7.1373027221003147, -0.0000000000002493, 22.6949581776302516, (-46.1374835353966546), (-0.0000000000004230), 35.1097108157676558});

	auto particle = P.get("D1232pp");
	particle->user_data("form_factor", Form_Factor::identity);

	auto s_values = weighted_space(1.1, particle->mass() - particle->width(), particle->mass() + particle->width(), 2.1, 17);

	auto diagram = create_diagram(FORMAT("{} s", particle->name()), s_channel, VMP,
	                              {Proton, Pi_Plus},
	                              {particle},
	                              {Pi_Plus, Proton}
	);

	Feynman_Process process({diagram});

	auto result = process.scattering_amplitude(s_values, {6, 8});

	for( auto const& row : result ){
		for( auto const& p : row ){
			std::cout << p.to_string('x') << "\n";
		}
		std::cout << "\n\n\n";
	}



//	s.save("s.poly");
//	c.save("c.poly");
//
//	Polynomial x, y;
//	x.load("s.poly");
//	y.load("c.poly");
//
//	std::cout << "s: " << s(1.3) << "\n";
//	std::cout << "x: " << x(1.3) << "\n";
//	std::cout << "c: " << c(0.4) << "\n";
//	std::cout << "y: " << y(0.4) << "\n";
//
//
//	func_t<1> s_func = s;
//	FPolynomial<1> p(s_func, c);


	/*
	double const start = cmd.as_double("start");
	double const end   = cmd.as_double("end");
	int    const steps = cmd.as_int("steps");


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
	*/
}