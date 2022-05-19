#ifndef EFFECTIVE_LAGRANGIAN_MODEL_HPP
#define EFFECTIVE_LAGRANGIAN_MODEL_HPP

#include <feynumeric/particle.hpp>
#include <feynumeric/particle_manager.hpp>
#include <feynumeric/vertex_manager.hpp>

extern Feynumeric::Vertex_Manager_Ptr VMP;
void init_vertices(Feynumeric::Particle_Manager const& P);

#endif // EFFECTIVE_LAGRANGIAN_MODEL_HPP