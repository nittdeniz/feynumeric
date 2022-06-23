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
				output_proc1 << sqrt_s;
				for( auto& elem : result ){
					output_proc1 << "\t" << elem;
				}
				output_proc1 << "\n";
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
			output_proc5 << sqrt_s;
			for( auto& [cost, row] : results ){
				output_proc5 << "\t" << cost;
				for( auto& elem : row ){
					output_proc5 << "\t" << elem;
				}
				output_proc5 << "\n";
			}
		};
	}
	stopwatch.stop();
	std::cout << "Finished. Time: " << stopwatch.time<std::chrono::milliseconds>()/1000. << "\n";
	std::string file_name = FORMAT("amplitude_{}", cmd.as_string("particle"));
	std::ofstream out(FORMAT("data/", file_name));
	out << "#decay to Npi\n";
	out << output_proc1.rdbuf();
	out << "#elastic scattering\n";
	out << output_proc5.rdbuf();
}

