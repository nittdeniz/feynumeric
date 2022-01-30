#include <feynumeric/dirac.hpp>
#include <feynumeric/units.hpp>
#include "particles.hpp"

using namespace Feynumeric::Units;
using Feynumeric::Particle;
using Feynumeric::Particle_Ptr;
using Feynumeric::Matrix;

Particle_Ptr Photon       = std::make_shared<Particle>("Photon", Particle::Type::Majorana, 0, 0, 1);

Particle_Ptr Electron     = std::make_shared<Particle>("Electron", Particle::Type::Particle, 0.5109989461_MeV, -1, 0.5);
Particle_Ptr Positron     = std::make_shared<Particle>("Positron", Particle::Type::AntiParticle, 0.5109989461_MeV, 1, 0.5);

Particle_Ptr Muon_Plus    = std::make_shared<Particle>("Muon_+", Particle::Type::AntiParticle, 105.6583745_MeV, 1, 0.5);
Particle_Ptr Muon_Minus   = std::make_shared<Particle>("Muon_-", Particle::Type::Particle, 105.6583745_MeV, -1, 0.5);

Particle_Ptr Proton       = std::make_shared<Particle>("Proton", Particle::Type::Particle, 938.270813_MeV, +1, 0.5);
Particle_Ptr Neutron      = std::make_shared<Particle>("Neutron", Particle::Type::Particle, 939.5654133_MeV, 0, 0.5);

Particle_Ptr Pi_Zero      = std::make_shared<Particle>("Pi_0", Particle::Type::Majorana, 134.9768_MeV, 0, 0);
Particle_Ptr Pi_Plus      = std::make_shared<Particle>("Pi_+", Particle::Type::AntiParticle, 139.57039_MeV, +1, 0);
Particle_Ptr Pi_Minus     = std::make_shared<Particle>("Pi_-", Particle::Type::Particle, 139.57039_MeV, -1, 0);

Particle_Ptr N1440p       = std::make_shared<Particle>("N1440_+", Particle::Type::Particle, 1.44_GeV, +1, 0.5);
Particle_Ptr N1440n       = std::make_shared<Particle>("N1440_0", Particle::Type::Particle, 1.44_GeV, 0, 0.5);


void init_particles()
{
    using namespace Feynumeric;
    Photon->feynman_virtual = [](Edge const* e){
        if( e == nullptr )
        {
            critical_error("Photon virtual edge is nullptr.");
        }
        return Projector(e->particle(), e->four_momentum(), e->get_lorentz_indices());
    };
    Photon->feynman_incoming = [](Edge const* e){return Matrix(1,1, 13);};
    Photon->feynman_outgoing = [](Edge const* e){return Matrix(any_cast<double>(Proton->user_data("coupling.ppgamma")));};
    Photon->user_data("coupling.ppgammma", 0.3);

    Electron->feynman_virtual = [](Edge const* e){return Propagator(e->particle(), e->four_momentum(), e->get_lorentz_indices());};
    Electron->feynman_incoming = [](Edge const* e){return u(e->four_momentum(), e->spin(), e->get_lorentz_indices());};
    Electron->feynman_outgoing = [](Edge const* e){return ubar(e->four_momentum(), e->spin(), e->get_lorentz_indices());};

//    Positron->feynman_virtual = [](){return Matrix(1,1,1);};
//    Positron->feynman_incoming = [](){return Matrix(1,4,1);};
//    Positron->feynman_outgoing = [](){return Matrix(4,1,1);};

    Muon_Minus->feynman_virtual = [](Edge const* e){return Propagator(e->particle(), e->four_momentum(), e->get_lorentz_indices());};
    Muon_Minus->feynman_incoming = [](Edge const* e){return u(e->four_momentum(), e->spin(), e->get_lorentz_indices());};
    Muon_Minus->feynman_outgoing = [](Edge const* e){return ubar(e->four_momentum(), e->spin(), e->get_lorentz_indices());};

//    Muon_Plus->feynman_virtual = [](){return Matrix(1,1,1);};
//    Muon_Plus->feynman_incoming = [](){return Matrix(1,4,1);};
//    Muon_Plus->feynman_outgoing = [](){return Matrix(4,1,1);};

//    Proton->feynman_virtual = [](){return Matrix(4, 4, {1,2,3,4, 5,6,7,8, 9,10,11,12, 13,14,15,16});};
//    Proton->feynman_incoming = [&](){return Matrix(u12(Matrix(4, 1, {1,2,3,4}), Angular_Momentum(0.5, 0.5)));};
//    Proton->feynman_outgoing = [](){return Matrix(ubar12(Matrix(1, 4, {1,2,3,4}), Angular_Momentum(0.5, 0.5)));};

//    Pi_Zero->feynman_virtual = [](){return Matrix(1,1,1);};
//    Pi_Zero->feynman_incoming = [](){return Matrix(1,1,1);};
//    Pi_Zero->feynman_outgoing = [](){return Matrix(1,1,1);};

//    N1440p->feynman_virtual = [](){return Matrix(4, 4, {1,2,3,4, 5,6,7,8, 9,10,11,12, 13,14,15,16});};
//    N1440p->feynman_incoming = [&](){return Matrix(u12(Matrix(4, 1, {1,2,3,4}), Angular_Momentum(0.5, 0.5)));};
//    N1440p->feynman_outgoing = [](){return Matrix(ubar12(Matrix(1, 4, {1,2,3,4}), Angular_Momentum(0.5, 0.5)));};
}