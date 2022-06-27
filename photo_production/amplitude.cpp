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

	particle->feynman_virtual = [](std::shared_ptr<Graph_Edge> e, Kinematics const& kin){ return Projector(e, kin); };
	particle->user_data("form_factor", identity);

	init_vertices(P, cmd.as_string("coupling_constants"));
	Timer stopwatch;
	stopwatch.start();

	std::vector<double> cos_theta_values;
	cos_theta_values.reserve(201);
	cos_theta_values.push_back(-0.999);
	for( std::size_t i = 1; i < 100; ++i ){
		cos_theta_values.push_back(-1 + 0.02 * i);
	}
	cos_theta_values.push_back(0.999);

	std::vector<std::vector<Point>> results1;
	std::vector<std::map<double, std::vector<Point>>> results5;
	results1.resize(particle->spin().n_states() * Proton->spin().n_states());
	results5.resize(2*Proton->spin().n_states());


	auto coupl_str = coupling_string(cmd.as_string("particle"), "N", "Pion");
	couplings.set(coupl_str, 1.);

	#pragma omp parallel for
	for( int i = 0; i < steps; ++i ){
		double sqrt_s = start + i * (end-start)/steps;
		if( cmd.as_string("particle") == "D1232" ){
			auto dummy = std::make_shared<Particle>(*particle);
			auto diagram_pi1 = create_diagram(FORMAT("decay {} to proton pi+", dummy->name()), Decay_1_to_2, VMP,
		                                  {dummy},
		                                  {},
                          {Proton, Pi_Plus}
			);
			Feynman_Process processes({diagram_pi1});
			dummy->mass(sqrt_s);
			auto result = processes.decay_M();
			#pragma omp critical
			{
				for( std::size_t j = 0; j < result.size(); ++j ){
					results1[j].push_back(Point{sqrt_s, result[j]});
				}
			};

		}
		if( cmd.as_string("particle") == "D1600" ){
			auto dummy = std::make_shared<Particle>(*particle);
			auto diagram_pi1 = create_diagram(FORMAT("decay {} to proton pi+", dummy->name()), Decay_1_to_2, VMP,
			                                  {dummy},
			                                  {},
			                                  {Proton, Pi_Plus}
			);
			Feynman_Process processes({diagram_pi1});
			dummy->mass(sqrt_s);
			auto result = processes.decay_M();
			#pragma omp critical
			{
				for( std::size_t j = 0; j < result.size(); ++j ){
					results1[j].push_back(Point{sqrt_s, result[j]});
				}
			};

		}
		std::map<double, std::vector<Complex>> results;
		for( auto cos_theta : cos_theta_values ){
			auto scattering = create_diagram(FORMAT("{} s", particle->name()), s_channel, VMP,
			                                 {Proton, Pi_Plus},
			                                 {particle},
			                                 {Pi_Plus, Proton}
			);
			Feynman_Process process({scattering});
			results[cos_theta] = process.M_costheta(sqrt_s, cos_theta);
		}
		#pragma omp critical
		{
			for( auto& [cost, row] : results ){
				for( std::size_t j = 0; j < row.size(); ++j ){
					results5[j][sqrt_s].push_back(Point{cost,row[j]});
				}
			}
		};
	}
	stopwatch.stop();
	std::cout << "Finished. Time: " << stopwatch.time<std::chrono::milliseconds>()/1000. << "\n";
	std::string file_name1 = FORMAT("data/amplitude_{}_{}.txt", cmd.as_string("particle"), "N_pi");
	std::string file_name5 = FORMAT("data/amplitude_{}_{}.txt", cmd.as_string("particle"), "scattering");
	std::ofstream out1(file_name1);

	out1 << "{";
	bool outer_first = true;
	for( auto& item : results1 ){
		if( !outer_first ){
			out1 << ",";
		}
		out1 << "{";
		bool inner_first = true;
		for( auto& point : item ){
			if( !inner_first ){
				out1 << ",";
			}
			out1 << FORMAT("{{{:f},{:f}+{:f}I}}", point.x, point.y.real(), point.y.imag());
			inner_first = false;
		}
		out1 << "}";
		outer_first = false;
	}
	out1 << "}\n";
	for( std::size_t i = 0; i <= 10; ++i ){
		for( auto& item : results1 ){
			Polynomial p(i);
			p.fit(item);
			out1 << p.to_string('x') << "\n";
		}
		out1 << "\n";
	}

	/*
	out1 << "{";
	bool outer_first = true;
	for( auto& item : results1 ){
		if( !outer_first ){
			out1 << ",";
		}
		out1 << "Fit[{";
		bool first = true;
		for( auto& [key, value] : item ){
			if( !first ){
				out1 << ",";
			}
			out1 << FORMAT("{{{:f},{:f}+I {:f}}}", key, value.real(), value.imag());
			first = false;
		}
		out1 << "},{1,s,s^2,s^3,s^4,s^5,s^6,s^7,s^8},s]";
		outer_first = false;
	}
	out1 << "}//Chop";

	std::ofstream out5(file_name5);

	out5 << "{";
	outer_first = true;
	for( auto& item : results5 ){
		if( !outer_first ){
			out5 << ",";
		}
		out5 << "{";
		bool first = true;
		for( auto& [sqrts, row] : item ){
			if( !first ){
				out5 << ",";
			}
			out5 << FORMAT("{{{:f},Fit[{{", sqrts);
			bool inner_first = true;
			for( auto& [cos, val] : row ){
				if( !inner_first ){
					out5 << ",";
				}
				out5 << FORMAT("{{{:f},{:f} + I {:f}}}", cos, val.real(), val.imag());
				inner_first = false;
			}
			first = false;
			out5 << "},{1,c,c^2,c^3,c^4,c^5,c^6,c^7,c^8},c]}";
		}
		outer_first = false;
		out5 << "}";
	}
	out5 << "}//Chop";
	 */
}

