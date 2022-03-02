/**
#ifndef PARTICLES_HPP
#define PARTICLES_HPP

#include <Feynumeric/particle.hpp>

// Mesons
extern Feynumeric::Particle_Ptr Pi_Zero;
extern Feynumeric::Particle_Ptr Pi_Minus;
extern Feynumeric::Particle_Ptr Pi_Plus;

// Gauge Bosons
extern Feynumeric::Particle_Ptr Photon;

// Leptons
extern Feynumeric::Particle_Ptr Electron;
extern Feynumeric::Particle_Ptr Positron;
extern Feynumeric::Particle_Ptr Muon_Plus;
extern Feynumeric::Particle_Ptr Muon_Minus;

// Baryons
extern Feynumeric::Particle_Ptr Proton;
extern Feynumeric::Particle_Ptr Neutron;

extern Feynumeric::Particle_Ptr N1440p;
extern Feynumeric::Particle_Ptr N1440n;
extern Feynumeric::Particle_Ptr N1520p;
extern Feynumeric::Particle_Ptr N1520p;
extern Feynumeric::Particle_Ptr N1535p;
extern Feynumeric::Particle_Ptr N1535n;

extern Feynumeric::Particle_Ptr D1232pp;
extern Feynumeric::Particle_Ptr D1232p;
extern Feynumeric::Particle_Ptr D1232n;
extern Feynumeric::Particle_Ptr D1232m;

void init_particles();

#endif // PARTICLES_HPP
 **/