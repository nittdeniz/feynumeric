#ifndef FEYNUMERIC_FEYNMAN_DIAGRAM_HPP
#define FEYNUMERIC_FEYNMAN_DIAGRAM_HPP

#include <vector>

#include "angular_momentum.hpp"
#include "feynman_graph.hpp"
#include "lorentz_index.hpp"
#include "vertex_manager.hpp"

namespace Feynumeric{
	class Feynman_Diagram{
	public:
		using Feynman_Rule = std::function<Matrix(Kinematics const& kin)>;
		using Momentum_Func = std::function<Four_Momentum(Particle_Ptr const&, Kinematics const&)>;
	private:
		Feynman_Graph _graph;
		Vertex_Manager_Ptr _VMP;
		std::vector<Lorentz_Index_Ptr> _lorentz_indices;
		std::vector<Angular_Momentum_Ptr> _spins;
		std::vector<Feynman_Rule> _amplitude;
		std::vector<Momentum_Func> _four_momenta;

		void trace_fermion_line(Feynman_Graph::Edge_Ptr const& ptr, Direction const& start_direction);
		void trace_fermion_line(Feynman_Graph::Edge_Ptr const& origin, Feynman_Graph::Vertex_Ptr const& ptr, Direction const& start_direction);

		void fix_momenta();

		void add_spin(Feynman_Graph::Edge_Ptr const& edge_ptr);
		void add_lorentz_indices(Feynman_Graph::Edge_Ptr const& edge_ptr);

		Momentum_Func generate_four_momentum(Direction const& direction, std::size_t pos) const;

		void iterate_indices();
		void iterate_spins();

	public:
		Feynman_Diagram(Topology const& topology, Vertex_Manager_Ptr const& VMP, std::vector<Particle_Ptr> const& incoming_particles, std::vector<Particle_Ptr> const& virtual_particles, std::vector<Particle_Ptr> const& outgoing_particles);
		void generate_amplitude();
		Four_Momentum four_momentum(std::size_t index, Particle_Ptr const& P, Kinematics const& kin);
		double dsigma_dcos(Kinematics const& kin);
		Vertex_Manager_Ptr Vertex_Manager();
	};
}

#endif // FEYNUMERIC_FEYNMAN_DIAGRAM_HPP