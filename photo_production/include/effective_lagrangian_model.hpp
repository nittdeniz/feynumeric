#ifndef PARTICLES_HPP
#define PARTICLES_HPP

#include <feynumeric/particle.hpp>
#include <feynumeric/vertex_manager.hpp>

// Mesons
extern Feynumeric::Particle_Ptr Pi_Zero;
extern Feynumeric::Particle_Ptr Pi_Minus;
extern Feynumeric::Particle_Ptr Pi_Plus;


// Baryons
extern Feynumeric::Particle_Ptr Proton;
extern Feynumeric::Particle_Ptr Neutron;
//
extern Feynumeric::Particle_Ptr N1440p;
extern Feynumeric::Particle_Ptr N1440n;
extern Feynumeric::Particle_Ptr N1520p;
extern Feynumeric::Particle_Ptr N1520p;
extern Feynumeric::Particle_Ptr N1535p;
extern Feynumeric::Particle_Ptr N1535n;
//
extern Feynumeric::Particle_Ptr D1232pp;
extern Feynumeric::Particle_Ptr D1232p;
extern Feynumeric::Particle_Ptr D1232n;
extern Feynumeric::Particle_Ptr D1232m;

//
extern Feynumeric::Vertex_Manager_Ptr VMP;


//
void init_particles();
void init_vertices();
//
#endif // PARTICLES_HPP