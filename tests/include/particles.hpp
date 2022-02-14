#ifndef PARTICLES_HPP
#define PARTICLES_HPP

#include "Feynumeric/particle.hpp"

namespace Particles
{
    using Feynumeric::Particle_Ptr;
    extern Particle_Ptr Proton
                      , Neutron
                      , Photon
                      , Pi_Zero
                      , Pi_Plus
                      , Pi_Minus
                      , Rho_zero
                      , Rho_plus
                      , Rho_minus
                      , N1440p
                      , N1440n
                      , N1520p
                      , N1520n
                      , N1535p
                      , N1535n
                      , D1232;


}

#endif // PARTICLES_HPP
