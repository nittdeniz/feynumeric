#ifndef FEYNUMERIC_QED_HPP
#define FEYNUMERIC_QED_HPP

#include "particle.hpp"
#include "vertex.hpp"
#include "vertex_manager.hpp"

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

		extern Feynumeric::Vertex_Manager_Ptr VMP;

		void init_particles();
		void init_vertices();
	}
}

#endif // FEYNUMERIC_QED_HPP