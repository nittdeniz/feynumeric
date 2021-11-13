#include "feyncalc/diagram.hpp"

#include <iostream>

namespace Feyncalc
{
    Diagram::Diagram(Graph const& graph, vector<Particle_Ptr> &&incoming_list, vector<Particle_Ptr> &&virtual_list,
                               vector<Particle_Ptr> &&outgoing_list)
                               : _graph(graph)
                               , _incoming_particles(std::move(incoming_list))
                               , _virtual_particles(std::move(virtual_list))
                               , _outgoing_particles(std::move(outgoing_list))
    {

    }

    void Diagram::try_generate_amplitude()
    {
        assert_diagram_validity();

        using std::cerr;
        using std::abort;

        // check incoming particles for fermions
        for( const auto& P : _incoming_particles )
        {
            if( P->is_fermion() || P->is_anti_fermion() )
            {

            }
        }

    }

    Complex Diagram::calculate_amplitude(const double sqrt_s, const double cos_theta) const
    {
        const double sin_theta = std::sqrt(1 - cos_theta * cos_theta);
        return Feyncalc::Complex(sqrt_s * sin_theta);
    }

    #ifdef CATCH2_TESTING_ENABLED
    bool Diagram::assert_diagram_validity() const
    #else
    void Diagram::assert_diagram_validity() const
    #endif
    {
        #ifdef CATCH2_TESTING_ENABLED
        return true;
        #endif
    }

}
