#include "feynman_process.hpp"
#include "feynman_diagram.hpp"
#include "four_vector.hpp"
#include "integrate.hpp"
#include "units.hpp"

#include <iomanip>
#include <format.hpp>
#include <iostream>

namespace Feynumeric
{
	Feynman_Process::Feynman_Process(std::initializer_list<Feynman_Diagram_Ptr> list)
			: _diagrams(list){
		validate_diagram_compatibility();
		for( auto& diagram : _diagrams ){
			diagram->generate_amplitude();
		}
	}

	void Feynman_Process::add_diagram(Feynman_Diagram_Ptr diagram){
		diagram->generate_amplitude();
		_diagrams.push_back(diagram);
		validate_diagram_compatibility();
	}

	void Feynman_Process::validate_diagram_compatibility() const{
		if( _diagrams.empty()){
			return;
		}

		auto const& incoming_particles = _diagrams[0]->incoming_particles();
		auto const& outgoing_particles = _diagrams[0]->outgoing_particles();
		for( auto const& diagram : _diagrams ){
			auto temp_i = diagram->incoming_particles();
			auto temp_o = diagram->outgoing_particles();
			if( temp_i != incoming_particles || temp_o != outgoing_particles ){
				critical_error("Diagrams are not compatible.");
			}
		}
	}

	void Feynman_Process::print_dsigma_dcos_table(std::ostream& out, double sqrt_s, std::size_t steps){
		std::vector<double> values(steps + 1);
		double const delta = 2. / steps;
		for( std::size_t i = 0; i < steps; ++i ){
			values[i] = -1 + i * delta;
		}
		values[steps] = 1.;
		print_dsigma_dcos_table(out, sqrt_s, std::move(values));
	}

	void Feynman_Process::print_dsigma_dcos_table(std::ostream& out, double sqrt_s, double delta){
		std::size_t const steps = 2. / delta;
		std::vector<double> values(steps + 1);
		for( std::size_t i = 0; i < steps; ++i ){
			values[i] = -1 + i * delta;
		}
		values[steps] = 1.;
		print_dsigma_dcos_table(out, sqrt_s, std::move(values));
	}

	void Feynman_Process::print_dsigma_dcos_table(std::ostream& out, double sqrt_s, std::vector<double>&& values){
		// table header
		out << "cos";
		for( auto const& diagram : _diagrams ){
			out << FORMAT("\t{} ({},{})", diagram->_name, diagram->_phase.real(), diagram->_phase.imag());
		}
		out << "\tsum\n";

		// data
		auto table = dsigma_dcos_table(sqrt_s, std::move(values));
		for( auto const&[cosine, results] : table ){
			out << cosine << "\t";
			for( auto const& value : results ){
				out << std::setw(10) << std::setprecision(10) << value << "\t";
			}
			out << "\n";
		}
	}

	std::map<double, std::vector<double>>
	Feynman_Process::dsigma_dcos_table(double sqrt_s, std::vector<double>&& values){
		using namespace Feynumeric::Units;

		Kinematics kin(sqrt_s, 2, 2);

		auto const& incoming = _diagrams[0]->incoming_particles();
		auto const& outgoing = _diagrams[0]->outgoing_particles();

		auto qin = momentum(sqrt_s, incoming[0]->mass(), incoming[1]->mass());
		auto qout = momentum(sqrt_s, outgoing[0]->mass(), outgoing[1]->mass());

		kin.incoming(0, four_momentum(qin, incoming[0]->mass(), 1));
		kin.incoming(1, four_momentum(-qin, incoming[1]->mass(), 1));

		for( auto& diagram : _diagrams ){
//			diagram->generate_amplitude();
			diagram->reset_spins();
			diagram->reset_indices();
		}

		std::size_t const N_spins = [&](){
			std::size_t n = 1;
			for( auto const& j : _diagrams[0]->_spins ){
				n *= j->n_states();
			}
			return n;
		}();

		std::size_t const N_polarisations = [&](){
			std::size_t n = 1;
			for( auto const& p : _diagrams[0]->_graph._incoming ){
				n *= p->spin()->n_states();
			}
			return n;
		}();

		double const phase_space_factor =
				1. / N_polarisations * 1. / ( 32 * M_PI * kin.sqrt_s() * kin.sqrt_s()) * 1._hbarc * 1._hbarc *
				fm_to_mub * qout / qin;

		std::map<double, std::vector<double>> result;

		for( std::size_t k = 0; k < values.size(); ++k ){
			auto const& cos_theta = values[k];
			kin.outgoing(0, four_momentum(qout, outgoing[0]->mass(), cos_theta));
			kin.outgoing(1, four_momentum(-qout, outgoing[1]->mass(), cos_theta));

			std::vector<double> Ms_squared(_diagrams.size() + 1);

			for( std::size_t i = 0; i < N_spins; ++i ){
				Complex M{0, 0};
				for( std::size_t j = 0; j < _diagrams.size(); ++j ){
					auto const& temp = _diagrams[j]->evaluate_amplitude(kin);
					M += temp;
					Ms_squared[j] += ( temp * std::conj(temp)).real();
					_diagrams[j]->iterate_spins();
				}
				Ms_squared[_diagrams.size()] += ( M * std::conj(M)).real();
			}

			for( auto& value : Ms_squared ){
				value *= phase_space_factor;
			}
			result[cos_theta] = Ms_squared;
		}
		return result;
	}

	std::map<double, std::vector<double>> Feynman_Process::dsigma_dcos_table(double sqrt_s, std::size_t steps){
		std::vector<double> values(steps + 1);
		double const delta = 2. / steps;
		for( std::size_t i = 0; i < steps; ++i ){
			values[i] = -1 + i * delta;
		}
		values[steps] = 1.;
		return dsigma_dcos_table(sqrt_s, std::move(values));
	}

	std::map<double, std::vector<double>> Feynman_Process::dsigma_dcos_table(double sqrt_s, double delta){
		std::size_t const steps = 2. / delta;
		std::vector<double> values(steps + 1);
		for( std::size_t i = 0; i < steps; ++i ){
			values[i] = -1 + i * delta;
		}
		values[steps] = 1.;
		return dsigma_dcos_table(sqrt_s, std::move(values));
	}

	double Feynman_Process::no_check_dsigma_dcos(double sqrt_s, double cos_theta){
		using namespace Feynumeric::Units;

		Kinematics kin(sqrt_s, 2, 2);

		auto const& incoming = _diagrams[0]->incoming_particles();
		auto const& outgoing = _diagrams[0]->outgoing_particles();

		auto qin = momentum(sqrt_s, incoming[0]->mass(), incoming[1]->mass());
		auto qout = momentum(sqrt_s, outgoing[0]->mass(), outgoing[1]->mass());

		kin.incoming(0, four_momentum(qin, incoming[0]->mass(), 1));
		kin.incoming(1, four_momentum(-qin, incoming[1]->mass(), 1));

		for( auto& diagram : _diagrams ){
//			diagram->generate_amplitude();
			diagram->reset_spins();
			diagram->reset_indices();
		}

		std::size_t const N_spins = [&](){
			std::size_t n = 1;
			for( auto const& j : _diagrams[0]->_spins ){
				n *= j->n_states();
			}
			return n;
		}();

		std::size_t const N_polarisations = [&](){
			std::size_t n = 1;
			for( auto const& p : _diagrams[0]->_graph._incoming ){
				n *= p->spin()->n_states();
			}
			return n;
		}();

		double const phase_space_factor =
				1. / N_polarisations * 1. / ( 32 * M_PI * kin.sqrt_s() * kin.sqrt_s()) * 1._hbarc * 1._hbarc *
				fm_to_mub * qout / qin;

		std::map<double, std::vector<double>> result;


			kin.outgoing(0, four_momentum(qout, outgoing[0]->mass(), cos_theta));
			kin.outgoing(1, four_momentum(-qout, outgoing[1]->mass(), cos_theta));

			std::vector<double> Ms_squared(_diagrams.size() + 1);

			for( std::size_t i = 0; i < N_spins; ++i ){
				Complex M{0, 0};
				for( std::size_t j = 0; j < _diagrams.size(); ++j ){
					auto const& temp = _diagrams[j]->evaluate_amplitude(kin);
					M += temp;
					Ms_squared[j] += ( temp * std::conj(temp)).real();
					_diagrams[j]->iterate_spins();
				}
				Ms_squared[_diagrams.size()] += ( M * std::conj(M)).real();
			}

			for( auto& value : Ms_squared ){
				value *= phase_space_factor;
			}
			result[cos_theta] = Ms_squared;

		auto check = dsigma_dcos_table(sqrt_s, 0.1);

		return result[cos_theta].back();
	}

	double Feynman_Process::decay_width(){
		if( _diagrams[0]->_graph._incoming.size() != 1 || _diagrams[0]->_graph._outgoing.size() != 2 ){
			critical_error("Only single particle decays into a two particle final state are implemented.");
		}
		using namespace Feynumeric::Units;
		validate_diagram_compatibility();

		auto const& incoming = _diagrams[0]->incoming_particles();
		auto const& outgoing = _diagrams[0]->outgoing_particles();

		Kinematics kin(incoming[0]->mass(), 1, 2);

		auto const q = momentum(incoming[0]->mass(), outgoing[0]->mass(), outgoing[1]->mass());

		kin.incoming(0, four_momentum(0, incoming[0]->mass()));
		kin.outgoing(0, four_momentum(q, outgoing[0]->mass()));
		kin.outgoing(1, four_momentum(-q, outgoing[1]->mass()));

		double Ms_squared{0.};

		for( auto& diagram : _diagrams ){
//			diagram->generate_amplitude();
			diagram->reset_spins();
			diagram->reset_indices();
		}

		std::size_t const N_spins = [&](){
			std::size_t n = 1;
			for( auto const& j : _diagrams[0]->_spins ){
				n *= j->n_states();
			}
			return n;
		}();

		std::size_t const N_polarisations = [&](){
			std::size_t n = 1;
			for( auto const& p : _diagrams[0]->_graph._incoming ){
				n *= p->spin()->n_states();
			}
			return n;
		}();

		for( std::size_t i = 0; i < N_spins; ++i ){
			Complex M{0, 0};
			for( auto& diagram : _diagrams ){
				auto const& temp = diagram->evaluate_amplitude(kin);
				M += temp;
				diagram->iterate_spins();
			}
			Ms_squared += ( M * std::conj(M)).real();
		}
		return Ms_squared *= 1. / N_polarisations * q / ( 8 * M_PI * incoming[0]->mass() * incoming[0]->mass());
	}

	void Feynman_Process::print_sigma_table(std::ostream& out, std::vector<double>&& values){
		using namespace Feynumeric::Units;
		using namespace std::placeholders;

		for( auto const& sqrt_s : values ){
			auto f = std::bind(&Feynman_Process::no_check_dsigma_dcos, this, sqrt_s, _1);
			out << FORMAT("{}\t{}\n", sqrt_s, integrate(f, -1., 1., 1.e-2));
		}
	}
}