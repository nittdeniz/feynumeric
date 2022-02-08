#include <iostream>


#include <feynumeric/feynman_diagram.hpp>
#include <feynumeric/process.hpp>
#include <feynumeric/topology.hpp>
#include <feynumeric/units.hpp>
#include <feynumeric/qed.hpp>


#include "particles.hpp"
#include "vertices.hpp"


int main()
{
    using namespace Feynumeric;
    using namespace Feynumeric::Units;

    init_particles();
	init_vertices();

    Topology Double_Wrench({
           {0, 2, Direction::INCOMING},
           {1,2, Direction::INCOMING},
           {2,3, Direction::VIRTUAL},
           {3, 4, Direction::OUTGOING},
           {3, 5, Direction::OUTGOING}
   });

    Topology X_Man({
		                   {0, 2, Direction::INCOMING},
		                   {1, 3, Direction::INCOMING},
		                   {2, 3, Direction::VIRTUAL},
		                   {2, 4, Direction::OUTGOING},
		                   {3,5, Direction::OUTGOING}
		    });

    Feynman_Diagram e_muon(X_Man, VMP,{Electron, Muon_Minus}, {Photon}, {Electron, Muon_Minus});

    Kinematics kin;
    kin.sqrt_s = 1.49_GeV;
    kin.cosines.resize(1);
    kin.momenta.push_back(0.3_GeV);
    kin.momenta.push_back(0.1_GeV);

	e_muon.generate_amplitude();

    std::cout << "result: " << e_muon.dsigma_dcos(kin);
    return EXIT_SUCCESS;
}

