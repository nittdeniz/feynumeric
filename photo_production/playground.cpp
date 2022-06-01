#include <feynumeric/feynumeric.hpp>
#include <feynumeric/qed.hpp>
#include "effective_lagrangian_model.hpp"
#include "form_factors.hpp"
#include <iostream>
#include <iomanip>

int main(int argc, char** argv){
    using namespace Feynumeric;
    using namespace Feynumeric::Units;

    //init_particles();


    if( argc != 2 )
    {
        critical_error("Program expects filename as a single and only argument.");
    }
    Particle_Manager P((std::string(argv[1])));
    init_vertices(P);
    auto const& Proton = P.get("proton");
    auto const& Neutron = P.get("neutron");
    //auto const& Photon = P.get("photon");
    auto const& Pi_Zero = P.get("pi0");
    auto const& Pi_Minus = P.get("pi-");
    auto const& Pi_Plus = P.get("pi+");


    P.get("N1440p")->user_data("form_factor", identity);
    P.get("N1440p")->user_data("gD1232N1440pi", 1.0);
    P.get("D1232pp")->user_data("gD1232N1440pi", 1.0);
    P.get("D1232p")->user_data("gD1232N1440pi", 1.0);
    P.get("D1232n")->user_data("gD1232N1440pi", 1.0);
    P.get("D1232m`")->user_data("gD1232N1440pi", 1.0);
    P.get("N1440n")->user_data("form_factor", identity);
    P.get("N1440n")->user_data("gD1232N1440pi", 1.0);
    P.get("D1232pp")->user_data("form_factor", identity);
    P.get("D1232p")->user_data("form_factor", identity);
    P.get("D1232n")->user_data("form_factor", identity);
    P.get("D1232m")->user_data("form_factor", identity);

    auto channel_decay_N1440p_D1232_pi_1 = create_diagram(FORMAT("decay N1440 to Delta pi 1"), Decay_1_to_M2_1, VMP,
                                                          {P.get("N1440p")},
                                                          {P.get("D1232pp")},
                                                          {Proton, Pi_Plus, Pi_Minus}
    );
    auto channel_decay_N1440p_D1232_pi_2 = create_diagram(FORMAT("decay N1440 to Delta pi 1"), Decay_1_to_M2_1_cross, VMP,
                                                          {P.get("N1440p")},
                                                          {P.get("D1232n")},
                                                          {Proton, Pi_Plus, Pi_Minus}
    );

    auto channel_decay_N1440p_D1232_pi_3 = create_diagram(FORMAT("decay N1440 to Delta pi 1"), Decay_1_to_M2_1, VMP,
                                                          {P.get("N1440p")},
                                                          {P.get("D1232p")},
                                                          {Proton, Pi_Zero, Pi_Zero}
    );
    auto channel_decay_N1440p_D1232_pi_4 = create_diagram(FORMAT("decay N1440 to Delta pi 1"), Decay_1_to_M2_1_cross, VMP,
                                                          {P.get("N1440p")},
                                                          {P.get("D1232p")},
                                                          {Proton, Pi_Zero, Pi_Zero}
    );

    auto channel_decay_N1440p_D1232_pi_5 = create_diagram(FORMAT("decay N1440 to Delta pi 1"), Decay_1_to_M2_1, VMP,
                                                          {P.get("N1440p")},
                                                          {P.get("D1232p")},
                                                          {Neutron, Pi_Plus, Pi_Zero}
    );
    auto channel_decay_N1440p_D1232_pi_6 = create_diagram(FORMAT("decay N1440 to Delta pi 1"), Decay_1_to_M2_1_cross, VMP,
                                                          {P.get("N1440p")},
                                                          {P.get("D1232n")},
                                                          {Neutron, Pi_Plus, Pi_Zero}
    );

    auto channel_decay_N1440n_D1232_pi_1 = create_diagram(FORMAT("decay N1440 to Delta pi 1"), Decay_1_to_M2_1, VMP,
                                                          {P.get("N1440n")},
                                                          {P.get("D1232n")},
                                                          {Proton, Pi_Minus, Pi_Zero}
    );
    auto channel_decay_N1440n_D1232_pi_2 = create_diagram(FORMAT("decay N1440 to Delta pi 1"), Decay_1_to_M2_1_cross, VMP,
                                                          {P.get("N1440n")},
                                                          {P.get("D1232p")},
                                                          {Proton, Pi_Minus, Pi_Zero}
    );

    auto channel_decay_N1440n_D1232_pi_3 = create_diagram(FORMAT("decay N1440 to Delta pi 1"), Decay_1_to_M2_1, VMP,
                                                          {P.get("N1440n")},
                                                          {P.get("D1232m")},
                                                          {Neutron, Pi_Minus, Pi_Plus}
    );
    auto channel_decay_N1440n_D1232_pi_4 = create_diagram(FORMAT("decay N1440 to Delta pi 1"), Decay_1_to_M2_1_cross, VMP,
                                                          {P.get("N1440n")},
                                                          {P.get("D1232p")},
                                                          {Neutron, Pi_Minus, Pi_Plus}
    );

    auto channel_decay_N1440n_D1232_pi_5 = create_diagram(FORMAT("decay N1440 to Delta pi 1"), Decay_1_to_M2_1, VMP,
                                                          {P.get("N1440n")},
                                                          {P.get("D1232n")},
                                                          {Neutron, Pi_Zero, Pi_Zero}
    );
    auto channel_decay_N1440n_D1232_pi_6 = create_diagram(FORMAT("decay N1440 to Delta pi 1"), Decay_1_to_M2_1_cross, VMP,
                                                          {P.get("N1440n")},
                                                          {P.get("D1232n")},
                                                          {Neutron, Pi_Zero, Pi_Zero}
    );

    Feynman_Process decay_Np1({channel_decay_N1440p_D1232_pi_1, channel_decay_N1440p_D1232_pi_2});
    Feynman_Process decay_Np2({channel_decay_N1440p_D1232_pi_3, channel_decay_N1440p_D1232_pi_4});
    Feynman_Process decay_Np3({channel_decay_N1440p_D1232_pi_5, channel_decay_N1440p_D1232_pi_6});

    Feynman_Process decay_Nn1({channel_decay_N1440n_D1232_pi_1, channel_decay_N1440n_D1232_pi_2});
    Feynman_Process decay_Nn2({channel_decay_N1440n_D1232_pi_3, channel_decay_N1440n_D1232_pi_4});
    Feynman_Process decay_Nn3({channel_decay_N1440n_D1232_pi_5, channel_decay_N1440n_D1232_pi_6});

    auto w1 = decay_Np1.decay_width();
    auto w2 = decay_Np2.decay_width();
    auto w3 = decay_Np3.decay_width();
    auto w4 = decay_Nn1.decay_width();
    auto w5 = decay_Nn2.decay_width();
    auto w6 = decay_Nn3.decay_width();

    double const literature_value = P.get("N1440")->width() * P.get("N1440")->user_data<double>("branching_N_pipi_D1232");
    std::cout << FORMAT("g: {}\n", P.get("D1232pp")->user_data<double>("gRNpi"));
    std::cout << FORMAT("g: {}\n", P.get("D1232pp")->user_data<double>("gD1232N1440pi"));
    std::cout << FORMAT("literature_value: {}\n", literature_value);
    std::cout << FORMAT("w1: {} w2: {} w3: {}\n", w1, w2, w3);
    std::cout << FORMAT("w4: {} w5: {} w6: {}\n", w4, w5, w6);

    std::cout << FORMAT("g(N1440+ -> D1232): " ) << std::setw(10) << std::setprecision(10)<< std::sqrt(literature_value / ( w1 + w2 + w3 )) << "\n";
    std::cout << FORMAT("g(N1440n -> D1232): " ) << std::setw(10) << std::setprecision(10)<< std::sqrt(literature_value / ( w4 + w5 + w6 )) << "\n";
}