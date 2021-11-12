#ifndef FEYNCALC_DIAGRAM
#define FEYNCALC_DIAGRAM

#include <vector>

#include "particle.hpp"
#include "topology.hpp"

namespace Feyncalc
{
    using std::vector;

    template<Topology topology>
    class Diagram
    {
    private:
        vector<Particle_Ptr> _incoming_particles;
        vector<Particle_Ptr> _virtual_particles;
        vector<Particle_Ptr> _outgoing_particles;

    public:
        Diagram(vector<Particle_Ptr>&& incoming_list, vector<Particle_Ptr>&& virtual_list, vector<Particle_Ptr>&& outgoing_list);
    };
}
#endif // FEYNCALC_DIAGRAM