#include "particles.hpp"
#include "units.hpp"
namespace Feynumeric
{
	using namespace Feynumeric::Units;
	[[maybe_unused]] Particle_Ptr Electron = std::make_shared<Particle>("e-", Particle::Type::Particle, 0.511_keV, -1._e, 0.5);
	[[maybe_unused]] Particle_Ptr Positron = std::make_shared<Particle>("e+", Particle::Type::AntiParticle, 0.511_keV, 1._e,0.5);
	[[maybe_unused]] Particle_Ptr Muon_Minus = std::make_shared<Particle>("mu-", Particle::Type::Particle, 0.511_keV, -1._e,.5);
	[[maybe_unused]] Particle_Ptr Muon_Plus = std::make_shared<Particle>("mu+", Particle::Type::AntiParticle, 0.511_keV, 1._e,0.5);

	[[maybe_unused]] Particle_Ptr Photon = std::make_shared<Particle>("mu+", Particle::Type::Majorana, 0, 0, 1);
}

