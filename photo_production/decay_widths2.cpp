#include <feynumeric/feynumeric.hpp>
#include "effective_lagrangian_model.hpp"
#include "form_factors.hpp"

int main(int argc, char** argv)
{
    using namespace Feynumeric;
    using namespace Feynumeric::Units;
    using namespace std::placeholders;

    Command_Line_Manager cmd(argc, argv);
    cmd.register_command("particle_file", true, "file with particle parameters");
    cmd.register_command("coupling_constants", true, "file with coupling constants");
    cmd.register_command("start_m", std::string("1.0"), "starting point");
    cmd.register_command("end_m", std::string("3.0"), "end value");
    cmd.register_command("start_s", std::string("1.0"), "starting point");
    cmd.register_command("end_s", std::string("3.0"), "end value");
    cmd.register_command("steps_m", std::string("200"), "steps");
    cmd.register_command("steps_s", std::string("200"), "steps");
    cmd.register_command("particle", true, "particle name");

    cmd.crash_on_missing_mandatory_command();

    Particle_Manager P(cmd.as_string("particle_file"));
    auto const& Proton = P.get("proton");
    auto const& Neutron = P.get("neutron");
    auto const& Pi_Zero = P.get("pi0");
    auto const& Pi_Minus = P.get("pi-");
    auto const& Pi_Plus = P.get("pi+");

    init_vertices(P, cmd.as_string("coupling_constants"));

    std::vector<std::string> resonances = {"Fictional12+_32", "Fictional12-_32", "Fictional32+_32", "Fictional32-_32", "Fictional52+_32", "Fictional52-_32", "Fictional72+_32", "Fictional72-_32"};

    for( auto const& p_str : resonances )
    {

        auto particle_str = FORMAT("{}p", p_str);
        auto particle = P.get(particle_str);

        particle->user_data("form_factor", Form_Factor::identity);

        auto diag1 = create_diagram(FORMAT("Decay {} -> p pi0", particle->name()), Decay_1_to_2, VMP,
                                    {particle}, {}, {Proton, Pi_Zero});
        auto diag2 = create_diagram(FORMAT("Decay {} -> n pi+", particle->name()), Decay_1_to_2, VMP,
                                    {particle}, {}, {Neutron, Pi_Plus});
        Feynman_Process proc1({diag1});
        Feynman_Process proc2({diag2});




        for( auto const &value: values )
        {
            data.emplace_back(value, proc1.decay_width(value) + proc2.decay_width(value));
        }
        Polynomial p(6);
        p.fit(data);
        std::cout << FORMAT("{}: {} / {}\n", particle->name(), p(particle->mass()).real(), proc1.decay_width(particle->mass())+ proc2.decay_width(particle->mass()));
        p.save(FORMAT("widths/{}_width_N_Pi.poly", p_str));
        std::cout << p.to_string('s') << "\n";
    }
}