#include <iostream>

#include <Feynumeric/diagram.hpp>
#include <Feynumeric/graph.hpp>
#include <Feynumeric/process.hpp>
#include <Feynumeric/topologies.hpp>
#include <Feynumeric/units.hpp>

#include "particles.hpp"

int main()
{
    using namespace Feynumeric;
    using namespace Feynumeric::Units;

    using std::cout;

    init_particles();

    Diagram s_channel_N1440(Topology::Double_Wrench, {Photon, Proton}, {N1440p}, {Pi_Zero, Proton});
//    Diagram<Topology::Double_Wrench> s_channel_N1520({Photon, Proton}, {N1520p}, {Pi_Zero, Proton});
//    Diagram<Topology::Double_Wrench> s_channel_N1535({Photon, Proton}, {N1535p}, {Pi_Zero, Proton});
//
//    Diagram u_channel_N1440(Topology::Scissors, {Photon, Proton}, {N1440p}, {Pi_Zero, Proton});
//    Diagram u_channel_N1520(Topology::Scissors, {Photon, Proton}, {N1520p}, {Pi_Zero, Proton});
//    Diagram u_channel_N1535(Topology::Scissors, {Photon, Proton}, {N1535p}, {Pi_Zero, Proton});

    Process photo_production;
    photo_production.add_diagrams({
        s_channel_N1440
        });

    auto result = photo_production.dsigma_dcos_theta(1.49_GeV, 1);
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

