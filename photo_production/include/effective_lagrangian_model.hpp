#ifndef PARTICLES_HPP
#define PARTICLES_HPP

#include <feynumeric/particle.hpp>
#include <feynumeric/vertex_manager.hpp>

// Mesons
extern Feynumeric::Particle_Ptr Pi_Zero;
extern Feynumeric::Particle_Ptr Pi_Minus;
extern Feynumeric::Particle_Ptr Pi_Plus;

extern Feynumeric::Particle_Ptr Rho_Zero;
extern Feynumeric::Particle_Ptr Rho_Minus;
extern Feynumeric::Particle_Ptr Rho_Plus;

// Baryons
extern Feynumeric::Particle_Ptr Proton;
extern Feynumeric::Particle_Ptr Neutron;
//
extern Feynumeric::Particle_Ptr N1440p;
extern Feynumeric::Particle_Ptr N1440n;
extern Feynumeric::Particle_Ptr N1520p;
extern Feynumeric::Particle_Ptr N1520n;
extern Feynumeric::Particle_Ptr N1535p;
extern Feynumeric::Particle_Ptr N1535n;
extern Feynumeric::Particle_Ptr N1650p;
extern Feynumeric::Particle_Ptr N1650n;
extern Feynumeric::Particle_Ptr N1675p;
extern Feynumeric::Particle_Ptr N1675n;
extern Feynumeric::Particle_Ptr N1680p;
extern Feynumeric::Particle_Ptr N1680n;
extern Feynumeric::Particle_Ptr N1700p;
extern Feynumeric::Particle_Ptr N1700n;
extern Feynumeric::Particle_Ptr N1710p;
extern Feynumeric::Particle_Ptr N1710n;
extern Feynumeric::Particle_Ptr N1720p;
extern Feynumeric::Particle_Ptr N1720p;
extern Feynumeric::Particle_Ptr N1875p;
extern Feynumeric::Particle_Ptr N1875n;
extern Feynumeric::Particle_Ptr N1880p;
extern Feynumeric::Particle_Ptr N1880n;
extern Feynumeric::Particle_Ptr N1895p;
extern Feynumeric::Particle_Ptr N1895n;
extern Feynumeric::Particle_Ptr N1900p;
extern Feynumeric::Particle_Ptr N1900n;
extern Feynumeric::Particle_Ptr N2060p;
extern Feynumeric::Particle_Ptr N2060n;
extern Feynumeric::Particle_Ptr N2100p;
extern Feynumeric::Particle_Ptr N2100n;
extern Feynumeric::Particle_Ptr N2120p;
extern Feynumeric::Particle_Ptr N2120n;
extern Feynumeric::Particle_Ptr N2190p;
extern Feynumeric::Particle_Ptr N2190n;
extern Feynumeric::Particle_Ptr N2220p;
extern Feynumeric::Particle_Ptr N2220n;
extern Feynumeric::Particle_Ptr N2250p;
extern Feynumeric::Particle_Ptr N2250n;
extern Feynumeric::Particle_Ptr N2600p;
extern Feynumeric::Particle_Ptr N2600n;



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