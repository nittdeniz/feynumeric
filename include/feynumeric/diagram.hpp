#ifndef Feynumeric_DIAGRAM
#define Feynumeric_DIAGRAM

#include <variant>
#include <vector>

#include "edge.hpp"
#include "function_types.hpp"
#include "graph.hpp"
#include "lorentz_index.hpp"
#include "matrix.hpp"
#include "particle.hpp"
#include "vertex.hpp"

namespace Feynumeric
{
    class Vertex_Manager;
    class Diagram : public std::enable_shared_from_this<Diagram>
    {
    private:
        std::shared_ptr<Vertex_Manager> _vertex_manager;
        Graph _graph;
        std::vector<Particle_Ptr> _incoming_particles;
        std::vector<Particle_Ptr> _virtual_particles;
        std::vector<Particle_Ptr> _outgoing_particles;
        std::vector<std::function<Matrix()>> _amplitude;

        std::vector<Lorentz_Index> _lorentz_indices;
        std::vector<Angular_Momentum> _angular_momenta;

        #ifdef CATCH2_TESTING_ENABLED
        bool assert_diagram_validity() const;
        #else
        void assert_diagram_validity() const;
        #endif

        Edge_Id _starting_edge_id;
        void trace_fermion_line(Edge_Id current_edge_id);

        void add_vertex(Edge_Id a, Edge_Id b);
        void add_vertex(Vertex_Id vertex_id);

        std::vector<Edge_Id> _remaining_edge_ids;
        std::map<Vertex_Id, std::vector<Edge_Id>> _remaining_vertices;

        void register_lorentz_indices();
        void register_angular_momenta();

        void fix_momenta();



    public:
        Diagram(std::shared_ptr<Vertex_Manager> const& vertex_manager, Graph const& graph, std::vector<Particle_Ptr>&& incoming_list, std::vector<Particle_Ptr>&& virtual_list, std::vector<Particle_Ptr>&& outgoing_list);
        Diagram(Diagram const& diagram);
        Diagram& operator=(Diagram const& diagram);

        void generate_amplitude();

        std::size_t n_total_external() const;

        Complex calculate_amplitude(const double sqrt_s, const double cos_theta) const;


    };
}
#endif // Feynumeric_DIAGRAM