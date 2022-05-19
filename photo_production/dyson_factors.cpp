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
    {/// N1440 -> Npi and N1440 -> Npipi
        auto dummy = std::make_shared<Particle>(*P.get("N1440p"));
        dummy->user_data("form_factor", identity);
        auto diagram_pi = create_diagram(FORMAT("decay {} to proton pi0", dummy->name()), Decay_1_to_2, VMP,
                                      {dummy},
                                      {},
                                      {Proton, Pi_Zero}
        );
        auto diagram_pipi_1 = create_diagram(FORMAT("decay {} to proton pi0", dummy->name()), Decay_1_to_M2_1, VMP,
                                        {dummy},
                                        {P.get("D1232pp")},
                                        {Proton, Pi_Plus, Pi_Minus}
        );
        auto diagram_pipi_2 = create_diagram(FORMAT("decay {} to proton pi0", dummy->name()), Decay_1_to_M2_1, VMP,
                                             {dummy},
                                             {P.get("D1232n")},
                                             {Proton, Pi_Plus, Pi_Minus}
        );
        diagram_pipi_2->cross_outgoing(1, 2);
        auto diagram_pipi_3 = create_diagram(FORMAT("decay {} to proton pi0", dummy->name()), Decay_1_to_1_M2, VMP,
                                             {dummy},
                                             {P.get("f0_500")},
                                             {Proton, Pi_Plus, Pi_Minus}
        );

        Feynman_Process decay_pi({diagram_pi});
        Feynman_Process decay_pipi({diagram_pipi_1, diagram_pipi_2, diagram_pipi_3});

        double Gamma0_1 = decay_pi.decay_width();
        std::map<double, double> dyson_factor;
        for( int i = 0; i <= steps; ++i ){
            double value = start + (end - start) / steps * i;
            dummy->mass(value);
            dyson_factor[value] = decay_pi.decay_width() / Gamma0_1;
        }
        write_map(dyson_factor, "data/dyson_factors/N1440.txt");
    }
}