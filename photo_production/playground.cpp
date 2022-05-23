#include <feynumeric/feynumeric.hpp>
#include <feynumeric/qed.hpp>
#include <iostream>

int main(){
    using namespace Feynumeric;
    using namespace Feynumeric::Units;
    using namespace Feynumeric::QED;

    init_particles();
    init_vertices();


    status("t_channel\n");

    Feynman_Diagram_Ptr diagram_t_channel = create_diagram("diagram_t_channel", t_channel, VMP,
                                                           {Electron, Electron},
                                                           {Photon},
                                                           {Electron, Electron});

    status("u_channel\n");

    Feynman_Diagram_Ptr diagram_u_channel = create_diagram("diagram_u_channel", u_channel, VMP,
                                                           {Electron, Electron},
                                                           {Photon},
                                                           {Electron, Electron});


    Feynman_Process e_scattering({diagram_t_channel, diagram_u_channel});
    e_scattering.conversion_factor(1._2mubarn);

    std::stringstream out;

    double const cos_theta = 0.2134;

    auto result = e_scattering.dsigma_dcos_table( 500._MeV, {{cos_theta}});
}