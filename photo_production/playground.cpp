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
    P.get("N1440p")->user_data("gD1232", 1.0);
    P.get("N1440n")->user_data("form_factor", identity);
    P.get("N1440n")->user_data("gD1232", 1.0);
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
                                                          {P.get("D1232n")},
                                                          {Proton, Pi_Minus, Pi_Plus}
    );
    auto channel_decay_N1440p_D1232_pi_4 = create_diagram(FORMAT("decay N1440 to Delta pi 1"), Decay_1_to_M2_1_cross, VMP,
                                                          {P.get("N1440p")},
                                                          {P.get("D1232pp")},
                                                          {Proton, Pi_Minus, Pi_Plus}
    );

    auto channel_decay_N1440n_D1232_pi_1 = create_diagram(FORMAT("decay N1440 to Delta pi 1"), Decay_1_to_M2_1, VMP,
                                                          {P.get("N1440n")},
                                                          {P.get("D1232m")},
                                                          {Neutron, Pi_Minus, Pi_Plus}
    );
    auto channel_decay_N1440n_D1232_pi_2 = create_diagram(FORMAT("decay N1440 to Delta pi 1"), Decay_1_to_M2_1_cross, VMP,
                                                          {P.get("N1440n")},
                                                          {P.get("D1232p")},
                                                          {Neutron, Pi_Minus, Pi_Plus}
    );
    auto channel_decay_N1440n_D1232_pi_3 = create_diagram(FORMAT("decay N1440 to Delta pi 1"), Decay_1_to_M2_1, VMP,
                                                          {P.get("N1440n")},
                                                          {P.get("D1232p")},
                                                          {Neutron, Pi_Plus, Pi_Minus}
    );
    auto channel_decay_N1440n_D1232_pi_4 = create_diagram(FORMAT("decay N1440 to Delta pi 1"), Decay_1_to_M2_1_cross, VMP,
                                                          {P.get("N1440n")},
                                                          {P.get("D1232m")},
                                                          {Neutron, Pi_Plus, Pi_Minus}
    );

    Feynman_Process a1({channel_decay_N1440p_D1232_pi_1});
    Feynman_Process a2({channel_decay_N1440p_D1232_pi_2});
    Feynman_Process a3({channel_decay_N1440p_D1232_pi_3});
    Feynman_Process a4({channel_decay_N1440p_D1232_pi_4});
    Feynman_Process b1({channel_decay_N1440n_D1232_pi_1});
    Feynman_Process b2({channel_decay_N1440n_D1232_pi_2});
    Feynman_Process b3({channel_decay_N1440n_D1232_pi_3});
    Feynman_Process b4({channel_decay_N1440n_D1232_pi_4});

    Feynman_Process decay_Np1({channel_decay_N1440p_D1232_pi_1, channel_decay_N1440p_D1232_pi_2});
    Feynman_Process decay_Nn2({channel_decay_N1440n_D1232_pi_3, channel_decay_N1440n_D1232_pi_4});

    std::cout << FORMAT("a1: {}, a2: {}\n", a1.decay_width(), a2.decay_width());
    std::cout << FORMAT("a3: {}, a4: {}\n", a3.decay_width(), a4.decay_width());
    std::cout << FORMAT("b1: {}, b2: {}\n", b1.decay_width(), b2.decay_width());
    std::cout << FORMAT("b3: {}, b4: {}\n", b3.decay_width(), b4.decay_width());

    auto w1 = decay_Np1.decay_width();
//    std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
//    std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
//    std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
//    std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
    auto w5 = decay_Nn2.decay_width();

    std::cout << FORMAT("w1: {}\nw5: {}\n", w1, w5);
}