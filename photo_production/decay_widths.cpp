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
    cmd.register_command("start", std::string("1.0"), "starting point");
    cmd.register_command("end", std::string("3.0"), "end value");
    cmd.register_command("steps", std::string("200"), "steps");
    cmd.register_command("form_factor", Form_Factor::CMD_FORM_FACTOR_NONE, FORMAT("which form factor to use ({}, {}, {}, {}, {}, {})", Form_Factor::CMD_FORM_FACTOR_NONE, Form_Factor::CMD_FORM_FACTOR_CASSING, Form_Factor::CMD_FORM_FACTOR_CUTKOSKY, Form_Factor::CMD_FORM_FACTOR_MANLEY, Form_Factor::CMD_FORM_FACTOR_MONIZ, Form_Factor::CMD_FORM_FACTOR_BREIT_WIGNER));
    cmd.register_command("particle", true, "particle name");

    cmd.crash_on_missing_mandatory_command();

    Particle_Manager P(cmd.as_string("particle_file"));
    auto const& Proton = P.get("proton");
    auto const& Neutron = P.get("neutron");
    auto const& Pi_Zero = P.get("pi0");
    auto const& Pi_Minus = P.get("pi-");
    auto const& Pi_Plus = P.get("pi+");

    init_vertices(P, cmd.as_string("coupling_constants"));

    std::vector<std::string> const nucleon_resonances = {"N1440", "N1520", "N1535", "N1650", "N1675", "N1680", "N1700", "N1710", "N1720", "N1875", "N1880", "N1895", "N1900"};
    std::vector<std::string> const delta_resonances = {"D1232", "D1600", "D1620", "D1700", "D1750", "D1900", "D1905", "D1910", "D1920", "D1930", "D1940", "D1950"};

    std::vector<std::string> resonances;
    resonances.insert(resonances.end(), nucleon_resonances.cbegin(), nucleon_resonances.cend());
    resonances.insert(resonances.end(), delta_resonances.cbegin(), delta_resonances.cend());

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
        const int N = 30;
        auto values = weighted_space(Proton->mass() + Pi_Plus->mass(), particle->mass() - particle->width() / 2., particle->mass() + particle->width() / 2.,
                                     2.5, N);
        std::vector<Point> data;
        data.reserve(N + 1);
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