#include "feyncalc/diagram.hpp"

namespace Feyncalc
{
    template<Topology topology>
    Diagram<topology>::Diagram(vector<Particle_Ptr> &&incoming_list, vector<Particle_Ptr> &&virtual_list,
                               vector<Particle_Ptr> &&outgoing_list)
                               : _incoming_particles(std::move(incoming_list))
                               , _virtual_particles(std::move(virtual_list))
                               , _outgoing_particles(std::move(outgoing_list))
    {

    }



    template class Diagram<Topology::Double_Wrench>;


}
