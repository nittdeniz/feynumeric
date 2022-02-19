#include "feynman_process.hpp"
#include "feynman_diagram.hpp"
#include "four_vector.hpp"
#include "units.hpp"

#include <iomanip>
#include <format.hpp>

namespace Feynumeric
{
	Feynman_Process::Feynman_Process(std::initializer_list<Feynman_Diagram_Ptr> list)
	: _diagrams(list)
	{

	}

	void Feynman_Process::add_diagram(Feynman_Diagram_Ptr diagram)
	{
		_diagrams.push_back(diagram);
	}

	void Feynman_Process::validate_diagram_compatibility() const
	{
		if( _diagrams.size() == 0 )
		{
			return;
		}

		auto const& incoming_particles = _diagrams[0]->incoming_particles();
		auto const& outgoing_particles = _diagrams[0]->outgoing_particles();
		for( auto const& diagram : _diagrams )
		{
			auto temp_i = diagram->incoming_particles();
			auto temp_o = diagram->outgoing_particles();
			if( temp_i != incoming_particles || temp_o != outgoing_particles )
			{
				critical_error("Diagrams are not compatible.");
			}
		}
	}

	void Feynman_Process::dsigma_dcos_table(std::ostream& out, double sqrt_s, std::size_t steps)
	{
		std::vector<double> values(steps+1);
		double const delta = 2./steps;
		for( std::size_t i = 0; i < steps; ++i )
		{
			values[i] = -1 + i * delta;
		}
		values[steps] = 1.;
		dsigma_dcos_table(out, sqrt_s, std::move(values));
	}

	void Feynman_Process::dsigma_dcos_table(std::ostream& out, double sqrt_s, double delta)
	{
		std::size_t const steps = 2./delta;
		std::vector<double> values(steps+1);
		for( std::size_t i = 0; i < steps; ++i )
		{
			values[i] = -1 + i * delta;
		}
		values[steps] = 1.;
		dsigma_dcos_table(out, sqrt_s, std::move(values));
	}

	void Feynman_Process::dsigma_dcos_table(std::ostream& out, double sqrt_s, std::vector<double>&& values)
	{
		using namespace Feynumeric::Units;
		validate_diagram_compatibility();

		Kinematics kin(sqrt_s, 2, 2);

		auto const& incoming = _diagrams[0]->incoming_particles();
		auto const& outgoing = _diagrams[0]->outgoing_particles();

		auto qin  = momentum(sqrt_s, incoming[0]->mass(), incoming[1]->mass());
		auto qout = momentum(sqrt_s, outgoing[0]->mass(), outgoing[1]->mass());

		kin.incoming(0, four_momentum(qin, incoming[0]->mass(), 1, 0));
		kin.incoming(1, four_momentum(-qin, incoming[1]->mass(), 1, 0));

		for( auto& diagram : _diagrams )
		{
			diagram->generate_amplitude();
			diagram->reset_spins();
		}

		std::size_t const N_spins = [&](){
			std::size_t n = 1;
			for( auto const& j : _diagrams[0]->_spins )
			{
				n *= j->n_states();
			}
			return n;
		}();

		std::size_t const N_polarisations = [&](){
			std::size_t n = 1;
			for( auto const& p : _diagrams[0]->_graph._incoming )
			{
				n *= p->spin()->n_states();
			}
			return n;
		}();

		double const phase_space_factor =
				1./N_polarisations * 1./(32 * M_PI * kin.sqrt_s() * kin.sqrt_s()) * 1._hbarc * 1._hbarc * fm_to_mub * qout/qin;

//				n_pol_N * n_pol_Gamma * 1. / (64 * M_PI * M_PI * kin.s()) * hbarc * hbarc * fm_to_mub * kin.qOut / kin.qIn;
		out << "cos";
		for( auto const& diagram : _diagrams )
		{
			out << FORMAT("\t{} ({},{})", diagram->_name, diagram->_phase.real(), diagram->_phase.imag());
		}
		out  << "\tsum\n";
		for( std::size_t k = 0; k < values.size(); ++k )
		{
			auto const& cos_theta = values[k];
			out << cos_theta << "\t";
			kin.outgoing(0, four_momentum(qout, outgoing[0]->mass(), cos_theta, 0));
			kin.outgoing(1, four_momentum(-qout, outgoing[1]->mass(), cos_theta, 0));

			std::vector<double> Ms_squared(_diagrams.size() + 1);

			for( std::size_t i = 0; i < N_spins; ++i ){
				Complex M{0, 0};
				for( std::size_t j = 0; j < _diagrams.size(); ++j ){
					auto const& temp = _diagrams[j]->evaluate_amplitude(kin);
					M += temp;
					Ms_squared[j] += (temp * std::conj(temp)).real();
					_diagrams[j]->iterate_spins();
				}
				Ms_squared[_diagrams.size()] += (M * std::conj(M)).real();
			}

			for( auto const& value : Ms_squared )
			{
				out << std::setw(10) << std::setprecision(10) << phase_space_factor * value << "\t";
			}
			out << "\n";
		}
	}
}