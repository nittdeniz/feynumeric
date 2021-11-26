#include <feyncalc/dirac.hpp>
#include <feyncalc/particle.hpp>
#include <feyncalc/units.hpp>


#include "particles.hpp"

using namespace Feyncalc::Units;
using Feyncalc::Particle;
using Feyncalc::Matrix;

Feyncalc::Particle_Ptr Photon   = std::make_shared<Feyncalc::Particle>("Photon", Particle::Type::Majorana, 0, 0, 1);

Feyncalc::Particle_Ptr Proton   = std::make_shared<Feyncalc::Particle>("Proton", Particle::Type::Particle, 938.270813_MeV, +1, 0.5);
Feyncalc::Particle_Ptr Neutron  = std::make_shared<Feyncalc::Particle>("Neutron", Particle::Type::Particle, 939.5654133_MeV, 0, 0.5);

Feyncalc::Particle_Ptr Pi_Zero  = std::make_shared<Feyncalc::Particle>("Pi_0", Particle::Type::Majorana, 134.9768_MeV, 0, 0);
Feyncalc::Particle_Ptr Pi_Plus  = std::make_shared<Feyncalc::Particle>("Pi_+", Particle::Type::AntiParticle, 139.57039_MeV, +1, 0);
Feyncalc::Particle_Ptr Pi_Minus = std::make_shared<Feyncalc::Particle>("Pi_-", Particle::Type::Particle, 139.57039_MeV, -1, 0);

Feyncalc::Particle_Ptr N1440p   = std::make_shared<Feyncalc::Particle>("N1440_+", Particle::Type::Particle, 1.44_GeV, +1, 0.5);
Feyncalc::Particle_Ptr N1440n   = std::make_shared<Feyncalc::Particle>("N1440_0", Particle::Type::Particle, 1.44_GeV, 0, 0.5);


void init_particles()
{
    using namespace Feyncalc;
    Photon->feynman_virtual = [](){return Matrix(1,1,1);};
    Photon->feynman_incoming = [](){return Matrix(1,1, 1);};
    Photon->feynman_outgoing = [&](){return Matrix(any_cast<double>(Proton->user_data("coupling.ppgamma")));};
    Photon->user_data("coupling.ppgammma", 0.3);

    Proton->feynman_virtual = [](){return Matrix(4, 4, {1,2,3,4, 5,6,7,8, 9,10,11,12, 13,14,15,16});};
    Proton->feynman_incoming = [&](){return Matrix(u12(Matrix(4, 1, {1,2,3,4}), Angular_Momentum(0.5, 0.5)));};
    Proton->feynman_outgoing = [](){return Matrix(ubar12(Matrix(1, 4, {1,2,3,4}), Angular_Momentum(0.5, 0.5)));};

    Pi_Zero->feynman_virtual = [](){return Matrix(1,1,1);};
    Pi_Zero->feynman_incoming = [](){return Matrix(1,1,1);};
    Pi_Zero->feynman_outgoing = [](){return Matrix(1,1,1);};

    N1440p->feynman_virtual = [](){return Matrix(4, 4, {1,2,3,4, 5,6,7,8, 9,10,11,12, 13,14,15,16});};
    N1440p->feynman_incoming = [&](){return Matrix(u12(Matrix(4, 1, {1,2,3,4}), Angular_Momentum(0.5, 0.5)));};
    N1440p->feynman_outgoing = [](){return Matrix(ubar12(Matrix(1, 4, {1,2,3,4}), Angular_Momentum(0.5, 0.5)));};
}