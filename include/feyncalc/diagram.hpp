#ifndef FEYNCALC_DIAGRAM
#define FEYNCALC_DIAGRAM

#include <vector>

#include "matrix.hpp"
#include "graph.hpp"
#include "particle.hpp"

namespace Feyncalc
{
    using std::vector;

    class Diagram
    {
    private:
        Graph _graph;
        vector<Particle_Ptr> _incoming_particles;
        vector<Particle_Ptr> _virtual_particles;
        vector<Particle_Ptr> _outgoing_particles;
        vector<function<Matrix()>> _amplitude;

        #ifdef CATCH2_TESTING_ENABLED
        bool assert_diagram_validity() const;
        #else
        void assert_diagram_validity() const;
        #endif

        void trace_fermion_line(vector<Edge>& remaining_edges, Edge const& starting_edge, Edge const& current_edge);

    public:
        Diagram(Graph const& graph, vector<Particle_Ptr>&& incoming_list, vector<Particle_Ptr>&& virtual_list, vector<Particle_Ptr>&& outgoing_list);
        Diagram(Diagram const& diagram);
        Diagram& operator=(Diagram const& diagram);

        void try_generate_amplitude();
        void generate_amplitude();
        Complex calculate_amplitude(const double sqrt_s, const double cos_theta) const;


    };
}
#endif // FEYNCALC_DIAGRAM