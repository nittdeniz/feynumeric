#ifndef PARTICLES_HPP
#define PARTICLES_HPP

#include <feynumeric/particle.hpp>


namespace Feynumeric
{
	namespace QED
	{
		// Gauge Bosons
		extern Feynumeric::Particle_Ptr Photon;
		// Leptons
		extern Feynumeric::Particle_Ptr Electron;
		extern Feynumeric::Particle_Ptr Positron;
		extern Feynumeric::Particle_Ptr Muon_Plus;
		extern Feynumeric::Particle_Ptr Muon_Minus;

		extern Feynumeric::Particle_Ptr Electron_Neutrino;
		extern Feynumeric::Particle_Ptr Electron_Anti_Neutrino;
		extern Feynumeric::Particle_Ptr Muon_Neutrino;
		extern Feynumeric::Particle_Ptr Muon_Anti_Neutrino;

		void init_particles();
	}
}


#endif // PARTICLES_HPP