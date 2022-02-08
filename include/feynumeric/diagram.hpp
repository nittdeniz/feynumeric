/*
#ifndef Feynumeric_DIAGRAM
#define Feynumeric_DIAGRAM

#include <memory>
#include <variant>
#include <vector>

#include "edge.hpp"
#include "function_types.hpp"
#include "graph.hpp"
#include "lorentz_index.hpp"
#include "matrix.hpp"
#include "momentum.hpp"
#include "particle.hpp"
#include "vertex.hpp"

namespace Feynumeric
{
    class Vertex_Manager;
    class Diagram : public std::enable_shared_from_this<Diagram>
    {
    private:

        Vertex_Manager* _vertex_manager;
        std::shared_ptr<Graph> _graph;
        std::vector<Particle_Ptr> _incoming_particles;
        std::vector<Particle_Ptr> _virtual_particles;
        std::vector<Particle_Ptr> _outgoing_particles;

        std::vector<Momentum_Func> _momenta;

        std::vector<std::function<Matrix(Kinematics const&)>> _amplitude;

        std::list<Lorentz_Index_Ptr> _lorentz_indices;
        std::list<Angular_Momentum_Ptr> _angular_momenta;


        #ifdef CATCH2_TESTING_ENABLED
        bool assert_diagram_validity() const;
        #else
        void assert_diagram_validity() const;
        #endif

        Edge* _starting_edge_ptr;
        void trace_fermion_line(Edge* current_edge);

        void add_vertex(Edge* a, Edge* b);
        void add_vertex(std::size_t vertex_id);

        std::vector<Edge*> _remaining_edges;
        std::map<std::size_t, std::vector<Edge*>> _remaining_vertices;

        void register_lorentz_indices();
        void register_angular_momenta();

        void iterate_spins();
        void iterate_indices();

        void fix_momenta();

        void generate_momentum_functions();

        void link_edges_to_this();


    public:
        Diagram(Vertex_Manager* vertex_manager, Graph const& graph, std::vector<Particle_Ptr>&& incoming_list, std::vector<Particle_Ptr>&& virtual_list, std::vector<Particle_Ptr>&& outgoing_list);
        Diagram(Diagram const& diagram);
        Diagram& operator=(Diagram const& diagram);

        void generate_amplitude();
		Momentum_Func four_momentum(Matrix const& M);

        std::size_t n_total_external() const;

        Complex calculate_amplitude(const double sqrt_s, const double cos_theta);


    };
}
#endif // Feynumeric_DIAGRAM
 */