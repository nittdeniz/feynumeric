#include <fstream>
#include <iostream>

#include <feynumeric/feynumeric.hpp>
#include "effective_lagrangian_model.hpp"
#include "form_factors.hpp"
#include "write_map.hpp"

int main(int argc, char** argv)
{
    using namespace Feynumeric;
    using namespace Feynumeric::Units;

    Command_Line_Manager cmd(argc, argv);
    cmd.register_command("particle_file", true, "file with particle parameters");
    cmd.register_command("start", std::string("1.0"), "starting point");
    cmd.register_command("end", std::string("3.0"), "end value");
    cmd.register_command("steps", std::string("200"), "steps");

    cmd.crash_on_missing_mandatory_command();

    Particle_Manager P(cmd.as_string("particle_file"));
    auto const& Proton = P["proton"];
    auto const& Neutron = P["neutron"];
    auto const& Photon = P["photon"];
    auto const& Pi_Zero = P["pi0"];
    auto const& Pi_Minus = P["pi-"];
    auto const& Pi_Plus = P["pi+"];

    init_vertices(P);

    double const start = cmd.as_double("start");
    double const end   = cmd.as_double("end");
    int    const steps = cmd.as_int("steps");
    {/// D1232 -> Npi
        auto dummy = std::make_shared<Particle>(*P.get("D1232pp"));
        dummy->user_data("form_factor", identity);
        auto decay_1 = create_diagram(FORMAT("decay {} to proton pi+", dummy->name()), Decay_1_to_2, VMP,
                                      {dummy},
                                      {},
                                      {Proton, Pi_Plus}
        );
        Feynman_Process decay({decay_1});
        double Gamma0 = decay.decay_width();
        std::map<double, double> dyson_factor;
        for( int i = 0; i <= steps; ++i ){
            double value = start + (end - start) / steps * i;
            dummy->mass(value);
            dyson_factor[value] = decay.decay_width()/Gamma0;
        }
        write_map(dyson_factor, "data/dyson_factors/D1232.txt");
    }
    {/// N1440 -> Npi
        auto dummy = std::make_shared<Particle>(*P.get("N1440p"));
        dummy->user_data("form_factor", identity);
        auto decay_1 = create_diagram(FORMAT("decay {} to proton pi0", dummy->name()), Decay_1_to_2, VMP,
                                      {dummy},
                                      {},
                                      {Proton, Pi_Zero}
        );
        Feynman_Process decay({decay_1});
        double Gamma0 = decay.decay_width();
        std::map<double, double> dyson_factor;
        for( int i = 0; i <= steps; ++i ){
            double value = start + (end - start) / steps * i;
            dummy->mass(value);
            dyson_factor[value] = decay.decay_width()/Gamma0;
        }
        write_map(dyson_factor, "data/dyson_factors/N1440.txt");
    }
}