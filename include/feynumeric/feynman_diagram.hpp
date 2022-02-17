#ifndef FEYNUMERIC_FEYNMAN_DIAGRAM_HPP
#define FEYNUMERIC_FEYNMAN_DIAGRAM_HPP

//#define DEBUG_AMPLITUDE 1

#include <initializer_list>
#include <memory>
#include <vector>

#include "angular_momentum.hpp"
#include "feynman_graph.hpp"
#include "lorentz_index.hpp"
#include "vertex_manager.hpp"

namespace Feynumeric{
	class Feynman_Diagram{
	public:
		using Feynman_Rule = std::function<Matrix(Kinematics const& kin)>;
		using Momentum_Func = std::function<Four_Vector(Particle_Ptr const&, Kinematics const&)>;
	private:
		Feynman_Graph _graph;
		Vertex_Manager_Ptr _VMP;
		std::vector<Lorentz_Index_Ptr> _lorentz_indices;
		std::vector<Angular_Momentum_Ptr> _spins;
		std::vector<Feynman_Rule> _amplitude;
		std::vector<Momentum_Func> _four_momenta;
		Complex _phase{1,0};

		void trace_fermion_line(Feynman_Graph::Edge_Ptr const& ptr, Direction const& start_direction, Feynman_Graph::Vertex_Direction vertex_direction);
		void trace_fermion_line(Feynman_Graph::Edge_Ptr const& origin, Feynman_Graph::Vertex_Ptr const& ptr, Direction const& start_direction);

		void fix_momenta();

		void add_spin(Feynman_Graph::Edge_Ptr const& edge_ptr);
		void add_lorentz_indices(Feynman_Graph::Edge_Ptr const& edge_ptr);

//		Momentum_Func generate_four_momentum(Direction const& direction, std::size_t pos) const;

		Four_Vector four_momentum();

		void iterate_indices();
		void iterate_spins();

	public:
		Feynman_Diagram(Topology const& topology, Vertex_Manager_Ptr const& VMP, std::initializer_list<Particle_Ptr> const& incoming_particles, std::initializer_list<Particle_Ptr> const& virtual_particles, std::initializer_list<Particle_Ptr> const& outgoing_particles);
		void generate_amplitude();
		Four_Vector four_momentum(std::size_t index, Particle_Ptr const& P, Kinematics const& kin);
		double dsigma_dcos(Kinematics& kin);
		Vertex_Manager_Ptr Vertex_Manager();

		std::vector<Particle_Ptr> incoming_particles() const;
		std::vector<Particle_Ptr> outgoing_particles() const;

		Complex evaluate_amplitude(Kinematics const& kin);

		void reset_spins();

		void phase(Complex phi);

		friend class Feynman_Process;
	};

	using Feynman_Diagram_Ptr = std::shared_ptr<Feynman_Diagram>;

	Feynman_Diagram_Ptr operator*(Complex phi, Feynman_Diagram_Ptr& p);
}

#endif // FEYNUMERIC_FEYNMAN_DIAGRAM_HPP