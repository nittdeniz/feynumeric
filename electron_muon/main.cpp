#include <iostream>

#include <feyncalc/diagram.hpp>
#include <feyncalc/process.hpp>
#include <feyncalc/topologies.hpp>
#include <feyncalc/units.hpp>

#include "particles.hpp"
#include "vertices.hpp"


int main()
{
    using namespace Feyncalc;
    using namespace Feyncalc::Units;

    using std::cout;

    init_particles();
    init_vertices();

    if( VMP == nullptr )
    {
        std::cerr << "nullptr\n";
        abort();
    }

    Diagram electron_muon_s(VMP,
                            Topology::X_Man,
                            {Photon, Muon_Minus},
                            {Photon},
                            {Electron, Muon_Minus});

    Process electron_muon_scattering;
    electron_muon_scattering.add_diagrams({
        electron_muon_s
        });

    auto result = electron_muon_scattering.dsigma_dcos_theta(1.49_GeV, 1);
    for( auto const& data : result )
    {
        cout  << "data: " << data << "\t";
    }
    cout << "\n";
//
//    table({"cos_theta", "dsigma", "+/-"}, [&photo_production](double cos_theta){return photo_production.dsigma_dcos_theta(1.49_GeV, cos_theta);}, {-1,1,0.1});
//
//    table({"E_cm", "E_lab", "sigma", "+/-"}, [&photo_production](double sqrt_s){return photo_production.sigma(sqrt_s);}, {1.1_GeV, 2.0_GeV, 0.1_GeV});

    return EXIT_SUCCESS;
}

