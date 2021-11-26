#ifndef PARTICLES_HPP
#define PARTICLES_HPP

#include <feyncalc/particle.hpp>

// Mesons
extern Feyncalc::Particle_Ptr Pi_Zero;
extern Feyncalc::Particle_Ptr Pi_Minus;
extern Feyncalc::Particle_Ptr Pi_Plus;

// Gauge Bosons
extern Feyncalc::Particle_Ptr Photon;

// Leptons
extern Feyncalc::Particle_Ptr Electron;
extern Feyncalc::Particle_Ptr Positron;
extern Feyncalc::Particle_Ptr Muon_Plus;
extern Feyncalc::Particle_Ptr Muon_Minus;

// Baryons
extern Feyncalc::Particle_Ptr Proton;
extern Feyncalc::Particle_Ptr Neutron;

extern Feyncalc::Particle_Ptr N1440p;
extern Feyncalc::Particle_Ptr N1440n;
extern Feyncalc::Particle_Ptr N1520p;
extern Feyncalc::Particle_Ptr N1520p;
extern Feyncalc::Particle_Ptr N1535p;
extern Feyncalc::Particle_Ptr N1535n;

extern Feyncalc::Particle_Ptr D1232pp;
extern Feyncalc::Particle_Ptr D1232p;
extern Feyncalc::Particle_Ptr D1232n;
extern Feyncalc::Particle_Ptr D1232m;

void init_particles();

#endif // PARTICLES_HPP