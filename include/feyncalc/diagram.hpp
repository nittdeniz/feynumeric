#ifndef FEYNCALC_DIAGRAM
#define FEYNCALC_DIAGRAM

#include <vector>

#include "particle.hpp"
#include "graph.hpp"

namespace Feyncalc
{
    using std::vector;

    class Diagram
    {
    private:
        const Graph _graph;
        vector<Particle_Ptr> _incoming_particles;
        vector<Particle_Ptr> _virtual_particles;
        vector<Particle_Ptr> _outgoing_particles;
        vector<function<Complex(Kinematics const&)>> _amplitude;

        #ifdef CATCH2_TESTING_ENABLED
        bool assert_diagram_validity() const;
        #else
        void assert_diagram_validity() const;
        #endif


    public:
        Diagram(Graph const& graph, vector<Particle_Ptr>&& incoming_list, vector<Particle_Ptr>&& virtual_list, vector<Particle_Ptr>&& outgoing_list);
        void try_generate_amplitude();
        Complex calculate_amplitude(const double sqrt_s, const double cos_theta) const;


    };
}
#endif // FEYNCALC_DIAGRAM