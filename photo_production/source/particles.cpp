#include <feyncalc/units.hpp>
#include <feyncalc/particle.hpp>


#include "particles.hpp"

using namespace Feyncalc::Units;
using Feyncalc::Particle;

Feyncalc::Particle_Ptr Photon   = std::make_shared<Feyncalc::Particle>("Photon", Particle::Type::Majorana, 0, 0, 1);

Feyncalc::Particle_Ptr Proton   = std::make_shared<Feyncalc::Particle>("Proton", Particle::Type::Particle, 938.270813_MeV, +1, 0.5);
Feyncalc::Particle_Ptr Neutron  = std::make_shared<Feyncalc::Particle>("Neutron", Particle::Type::Particle, 939.5654133_MeV, 0, 0.5);

Feyncalc::Particle_Ptr Pi_Zero  = std::make_shared<Feyncalc::Particle>("Pi_0", Particle::Type::Majorana, 134.9768_MeV, 0, 0);
Feyncalc::Particle_Ptr Pi_Plus  = std::make_shared<Feyncalc::Particle>("Pi_+", Particle::Type::AntiParticle, 139.57039_MeV, +1, 0);
Feyncalc::Particle_Ptr Pi_Minus = std::make_shared<Feyncalc::Particle>("Pi_-", Particle::Type::Particle, 139.57039_MeV, -1, 0);

Feyncalc::Particle_Ptr N1440p   = std::make_shared<Feyncalc::Particle>("N1440_+", Particle::Type::Particle, 1.44_GeV, +1, 0.5);
Feyncalc::Particle_Ptr N1440n   = std::make_shared<Feyncalc::Particle>("N1440_0", Particle::Type::Particle, 1.44_GeV, 0, 0.5);

