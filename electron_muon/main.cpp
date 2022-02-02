#include <iostream>

#include <Feynumeric/diagram.hpp>
#include <Feynumeric/process.hpp>
#include <Feynumeric/topologies.hpp>
#include <Feynumeric/units.hpp>

#include "particles.hpp"
#include "vertices.hpp"


int main()
{
    using namespace Feynumeric;
    using namespace Feynumeric::Units;
    using std::cout;


    init_particles();
    init_vertices();

    Process electron_muon_scattering;

    Diagram electron_muon_s(&VM,
                            Topology::X_Man,
                            {Electron, Muon_Minus},
                            {Photon},
                            {Electron, Muon_Minus});



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

