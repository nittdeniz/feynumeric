#include "effective_lagrangian_model.hpp"
#include "form_factors.hpp"

#include <feynumeric/feynumeric.hpp>

#include <iostream>

int main(int argc, char** argv){
	using namespace Feynumeric;

	Command_Line_Manager cmd(argc, argv);
	cmd.register_command("particle_file", true, "file with particle parameters");
	cmd.register_command("coupling_constants", true, "file with coupling constants");
	cmd.register_command("values", true);
    cmd.register_command("topology", std::string("s-channel"));
	cmd.register_command("particle", true);
    cmd.register_command("in1", true);
    cmd.register_command("in2", true);
    cmd.register_command("out1", true);
    cmd.register_command("out2", true);
    cmd.register_command("sqrt_s", true);
	cmd.crash_on_missing_mandatory_command();

	Particle_Manager P(cmd.as_string("particle_file"));
	init_vertices(P, cmd.as_string("coupling_constants"));


	auto in1 = P[cmd.as_string("in1")];
    auto in2 = P[cmd.as_string("in2")];
    auto out1 = P[cmd.as_string("out1")];
    auto out2 = P[cmd.as_string("out2")];
    auto particle = P[cmd.as_string("particle")];

    particle->user_data("form_factor", Form_Factor::identity);
    particle->feynman_virtual = [](std::shared_ptr<Graph_Edge> e, Kinematics const& kin){ return Projector(e, kin); };

    couplings.set(coupling_string(cmd.as_string("in1"), cmd.as_string("in2"), cmd.as_string("particle")), 1.);
    couplings.set(coupling_string(cmd.as_string("in1"), cmd.as_string("out1"), cmd.as_string("particle")), 1.);
    couplings.set(coupling_string(cmd.as_string("in1"), cmd.as_string("out2"), cmd.as_string("particle")), 1.);
    couplings.set(coupling_string(cmd.as_string("in2"), cmd.as_string("out1"), cmd.as_string("particle")), 1.);
    couplings.set(coupling_string(cmd.as_string("in2"), cmd.as_string("out2"), cmd.as_string("particle")), 1.);
    couplings.set(coupling_string(cmd.as_string("out1"), cmd.as_string("out2"), cmd.as_string("particle")), 1.);

	std::vector<double> cos_values;
    {
        std::stringstream stream(cmd.as_string("values"));
        stream >> std::skipws;
        std::string buffer;
        while( std::getline(stream, buffer, ',')){
            try {
                double temp = std::stod(buffer);
                cos_values.push_back(temp);
            }
            catch( std::exception& e) {
                warning(std::string(e.what()));
            }
        }
    }

    Topology const* topology;
    if( cmd.as_string("topology") == "s-channel" ) topology = &s_channel;
    else if( cmd.as_string("topology") == "t-channel" ) topology = &t_channel;
    else critical_error("unknown or unsupported topology");

    double sqrt_s = cmd.as_double("sqrt_s");


    auto diagram = create_diagram("diagram", *topology, VMP, {in1, in2}, {particle}, {out1, out2});
    Feynman_Process process({diagram});

    class NullBuffer : public std::streambuf{
        public: int overflow(int c){return c;}
    } null_buffer;
    auto null_stream = std::ostream(&null_buffer);
    process.print_dsigma_dcos_table_trace(std::cout, null_stream, sqrt_s, std::move(cos_values));
    return EXIT_SUCCESS;
}