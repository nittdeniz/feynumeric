#include <omp.h>

#include <iomanip>
#include <iostream>

#include "dirac.hpp"
#include "feynman_process.hpp"
#include "feynman_diagram.hpp"
#include "format.hpp"
#include "four_vector.hpp"
#include "graph_edge.hpp"
#include "graph_vertex.hpp"
#include "integrate.hpp"
#include "particle.hpp"
#include "phase_space.hpp"
#include "utility.hpp"

void nothing(){}

namespace Feynumeric
{
	Feynman_Process::Feynman_Process(std::initializer_list<Feynman_Diagram_Ptr> list)
			: _diagrams(list){
		validate_diagram_compatibility();
		for( auto& diagram : _diagrams ){
			diagram->generate_amplitude();
		}
        _n_spins = n_spins();
		_n_polarisations = n_polarisations();
	}

	Feynman_Process::Feynman_Process(std::vector<Feynman_Diagram_Ptr> list)
			: _diagrams(std::move(list)){
		validate_diagram_compatibility();
		for( auto& diagram : _diagrams ){
			diagram->generate_amplitude();
		}
        _n_spins = n_spins();
        _n_polarisations = n_polarisations();
	}

    Feynman_Process::Feynman_Process(const Feynman_Process &other)
    : _conversion_factor(other._conversion_factor)
    , _n_spins(other._n_spins)
    , _n_polarisations(other._n_polarisations)
    {
        _diagrams.reserve(other._diagrams.size());
        for( auto const& diagram : other._diagrams ){
            _diagrams.push_back(std::make_shared<Feynman_Diagram>(*diagram));
            _diagrams.back()->generate_amplitude();
        }
    }

	Feynman_Process& Feynman_Process::operator=(Feynman_Process const& other){
		_conversion_factor = other._conversion_factor;
		_n_spins = other._n_spins;
		_n_polarisations = other._n_polarisations;
		_diagrams.reserve(other._diagrams.size());
		for( auto const& diagram : other._diagrams ){
			_diagrams.push_back(std::make_shared<Feynman_Diagram>(*diagram));
			_diagrams.back()->generate_amplitude();
		}
		return *this;
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
			out << FORMAT("\t{}", diagram->_name);
		}
		out << "\tsum\n";

		// data
		auto table = dsigma_dcos_table(sqrt_s, std::move(values));
		for( auto const&[cosine, results] : table ){
			out << cosine << "\t";
			for( auto const& value : results ){
				out << std::setw(10) << std::fixed << std::setprecision(10) << value << "\t";
			}
			out << "\n";
		}
	}

    std::size_t Feynman_Process::n_spins()
    {
        std::size_t n = 1;
        for( auto const& j : _diagrams[0]->_spins ){
            n *= j->n_states();
        }
        return n;
    }

    std::size_t Feynman_Process::n_polarisations()
    {
        std::size_t n = 1;
        for( auto const& p : _diagrams[0]->_graph._incoming ){
            n *= p->spin()->n_states();
        }
        return n;
    }

    void Feynman_Process::conversion_factor(long double x){
		_conversion_factor = static_cast<double>(x);
	}

	std::map<double, std::vector<double>>
	Feynman_Process::dsigma_dcos_table(double sqrt_s, std::vector<double>&& values){
		using namespace Feynumeric::Units;

        for( auto& diagram : _diagrams ){
            diagram->reset_spins();
            diagram->reset_indices();
        }

		Kinematics kin(sqrt_s, 2, 2);

		auto const& incoming = _diagrams[0]->incoming_particles();
		auto const& outgoing = _diagrams[0]->outgoing_particles();

		auto qin = momentum(sqrt_s, incoming[0]->mass(), incoming[1]->mass());
		auto qout = momentum(sqrt_s, outgoing[0]->mass(), outgoing[1]->mass());

		kin.incoming(0, four_momentum(qin, incoming[0]->mass(), 1));
		kin.incoming(1, four_momentum(-qin, incoming[1]->mass(), 1));

		std::map<double, std::vector<double>> result;

        double const phase_space_factor = phase_space2(_n_polarisations, kin.sqrt_s(), qout, qin);

		for( std::size_t k = 0; k < values.size(); ++k ){
			auto const& cos_theta = values[k];
			kin.outgoing(0, four_momentum(qout, outgoing[0]->mass(), cos_theta));
			kin.outgoing(1, four_momentum(-qout, outgoing[1]->mass(), cos_theta));

			std::vector<double> Ms_squared(_diagrams.size() + 1);

			for( std::size_t i = 0; i < _n_spins; ++i ){
				Complex M{0, 0};
				for( std::size_t j = 0; j < _diagrams.size(); ++j ){
					auto const& temp = _diagrams[j]->evaluate_amplitude(kin);
					M += temp;
					Ms_squared[j] += (temp * std::conj(temp)).real();
					_diagrams[j]->iterate_spins();
				}
				auto M2 = M * std::conj(M);
				Ms_squared[_diagrams.size()] += (M2).real();
			}
			for( auto& value : Ms_squared ){
				value *= phase_space_factor * _conversion_factor;
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

	std::vector<Complex> Feynman_Process::M_costheta(double sqrt_s, double cos_theta){
		using namespace Feynumeric::Units;

		Kinematics kin(sqrt_s, 2, 2);

		auto const& incoming = _diagrams[0]->incoming_particles();
		auto const& outgoing = _diagrams[0]->outgoing_particles();

		auto qin = momentum(sqrt_s, incoming[0]->mass(), incoming[1]->mass());
		auto qout = momentum(sqrt_s, outgoing[0]->mass(), outgoing[1]->mass());

		kin.incoming(0, four_momentum(qin, incoming[0]->mass(), 1));
		kin.incoming(1, four_momentum(-qin, incoming[1]->mass(), 1));

		for( auto& diagram : _diagrams ){
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

		kin.outgoing(0, four_momentum(qout, outgoing[0]->mass(), cos_theta));
		kin.outgoing(1, four_momentum(-qout, outgoing[1]->mass(), cos_theta));

		std::vector<Complex> Ms(N_spins);

		for( std::size_t i = 0; i < N_spins; ++i ){
			Complex M{0, 0};
			for( std::size_t j = 0; j < _diagrams.size(); ++j ){
				M += _diagrams[j]->evaluate_amplitude(kin);
				_diagrams[j]->iterate_spins();
			}
			Ms[i] = M;
		}
		return Ms;
	}

	double Feynman_Process::no_check_dsigma_dcos(double sqrt_s, double cos_theta){
		using namespace Feynumeric::Units;

		Kinematics kin(sqrt_s, 2, 2);
        kin.angle(0, cos_theta);

		auto const& incoming = _diagrams[0]->incoming_particles();
		auto const& outgoing = _diagrams[0]->outgoing_particles();

		auto qin = momentum(sqrt_s, incoming[0]->mass(), incoming[1]->mass());
		auto qout = momentum(sqrt_s, outgoing[0]->mass(), outgoing[1]->mass());

		kin.incoming(0, four_momentum(qin, incoming[0]->mass(), 1));
		kin.incoming(1, four_momentum(-qin, incoming[1]->mass(), 1));


		for( auto& diagram : _diagrams ){
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

		double const phase_space_factor = phase_space2(N_polarisations, kin.sqrt_s(), qout, qin);

//		std::cout << FORMAT("s: {}\t phi: {}\n", kin.sqrt_s(), phase_space_factor);

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
//            std::cout << "s: " << sqrt_s << " cos: " << cos_theta << " i: " << i << " M: " << M << "\n";
            auto M2 = ( M * std::conj(M)).real();
            Ms_squared[_diagrams.size()] += M2;
        }

        for( auto& value : Ms_squared ){
//        	std::cout << FORMAT("sqrts: {} costheta: {} value: {} ps: {} cv: {}\n", sqrt_s, cos_theta, value, phase_space_factor, _conversion_factor);
            value *= phase_space_factor * _conversion_factor;
        }
        result[cos_theta] = Ms_squared;

		//auto check = dsigma_dcos_table(sqrt_s, 0.1);

		return result[cos_theta].back();
	}

	std::pair<std::vector<Polynomial>, std::vector<Polynomial>> Feynman_Process::M(double from, double to, double norm, std::size_t order_cos, std::size_t order_sqrts){
		std::vector<Polynomial> cos_theta_polynomial = M_costheta_polynomial(norm, order_cos);

		auto s_values = lin_space(from, to, 2*order_sqrts);


		std::size_t const N_spins = [&](){
			std::size_t n = 1;
			for( auto const& j : _diagrams[0]->_spins ){
				n *= j->n_states();
			}
			return n;
		}();

		auto const cos_theta = 0.;

		std::vector<std::vector<Point>> point_list(N_spins, std::vector<Point>(s_values.size()));
		for( std::size_t k = 0; k < s_values.size(); ++k ){
			auto const& sqrt_s = s_values[k];

			Kinematics kin(sqrt_s, 2, 2);

			auto const& incoming = _diagrams[0]->incoming_particles();
			auto const& outgoing = _diagrams[0]->outgoing_particles();

			auto qin = momentum(sqrt_s, incoming[0]->mass(), incoming[1]->mass());
			auto qout = momentum(sqrt_s, outgoing[0]->mass(), outgoing[1]->mass());

			kin.incoming(0, four_momentum(qin, incoming[0]->mass(), 1));
			kin.incoming(1, four_momentum(-qin, incoming[1]->mass(), 1));

			for( auto& diagram : _diagrams ){
				diagram->reset_spins();
				diagram->reset_indices();
			}

			kin.outgoing(0, four_momentum(qout, outgoing[0]->mass(), cos_theta));
			kin.outgoing(1, four_momentum(-qout, outgoing[1]->mass(), cos_theta));

			std::vector<Complex> Ms(N_spins);

			for( std::size_t i = 0; i < N_spins; ++i ){
				Complex M{0, 0};
				for( auto& diagram : _diagrams ){
					M += diagram->evaluate_amplitude(kin);
					diagram->iterate_spins();
				}
				point_list[i][k] = {sqrt_s, M};
			}
		}

		std::vector<Polynomial> s_polynomials;
		s_polynomials.reserve(point_list.size());
		for( auto points : point_list ){
			Polynomial p(order_sqrts);
			p.fit(points);
//			std::cout << "evaluate_at: " << evaluate_cos_at << "\n";
			auto rescale = p(norm);
//			std::cout << "rescale: " << rescale << "\n";
			if( rescale == Complex(0., 0.) ){
				critical_error("rescale factor is zero. Try a different order for sqrt_s polynomial.");
			}
			std::cout << "rescale: " << rescale << "\n\n";
			s_polynomials.push_back(p/rescale); // rescale
//			std::cout << "{";
//			for( auto& pp : points ){
//				std::cout << "{" << pp.x << ", " << pp.y << "},";
//			}
//			std::cout << "}\n";
//			std::cout << p.to_string('c') << "\n\n";
		}
		return std::make_pair(cos_theta_polynomial, s_polynomials);
	}

	std::vector<Polynomial> Feynman_Process::M_costheta_polynomial(double sqrt_s, std::size_t order){
		using namespace Feynumeric::Units;

		auto cos_values = lin_space(-0.999,0.999, 2*order);

		std::size_t const N_spins = [&](){
			std::size_t n = 1;
			for( auto const& j : _diagrams[0]->_spins ){
				n *= j->n_states();
			}
			return n;
		}();

		std::vector<std::vector<Point>> point_list(N_spins, std::vector<Point>(cos_values.size()));

		for( std::size_t k = 0; k < cos_values.size(); ++k ){
			auto const& cos_theta = cos_values[k];
			Kinematics kin(sqrt_s, 2, 2);

			auto const& incoming = _diagrams[0]->incoming_particles();
			auto const& outgoing = _diagrams[0]->outgoing_particles();

			auto qin = momentum(sqrt_s, incoming[0]->mass(), incoming[1]->mass());
			auto qout = momentum(sqrt_s, outgoing[0]->mass(), outgoing[1]->mass());

			kin.incoming(0, four_momentum(qin, incoming[0]->mass(), 1));
			kin.incoming(1, four_momentum(-qin, incoming[1]->mass(), 1));

			for( auto& diagram : _diagrams ){
				diagram->reset_spins();
				diagram->reset_indices();
			}

			kin.outgoing(0, four_momentum(qout, outgoing[0]->mass(), cos_theta));
			kin.outgoing(1, four_momentum(-qout, outgoing[1]->mass(), cos_theta));

			std::vector<Complex> Ms(N_spins);

			for( std::size_t i = 0; i < N_spins; ++i ){
				Complex M{0, 0};
				for( auto& diagram : _diagrams ){
					M += diagram->evaluate_amplitude(kin);
					diagram->iterate_spins();
				}
				point_list[i][k] = {cos_theta, M};
			}
		}
		std::vector<Polynomial> polynomials;
		polynomials.reserve(point_list.size());
		for( auto points : point_list ){
			Polynomial p(order);
			p.fit(points);
			polynomials.push_back(p);
//			std::cout << "{";
//			for( auto& pp : points ){
//				std::cout << "{" << pp.x << ", " << pp.y << "},";
//			}
//			std::cout << "}\n";
//			std::cout << p.to_string('c') << "\n\n";
		}
		return polynomials;
	}

	std::vector<Polynomial> Feynman_Process::decay_M_polynomial(Particle_Ptr dummy, double from, double to, std::size_t order){
		validate_diagram_compatibility();

		std::size_t const N_spins = [&](){
			std::size_t n = 1;
			for( auto const& j : _diagrams[0]->_spins ){
				n *= j->n_states();
			}
			return n;
		}();

		auto x2 = dummy->mass() - dummy->width();
		auto x3 = dummy->mass() + dummy->width();

		std::vector<double> values1 = lin_space(from, x2, order/2-1);
		std::vector<double> values2 = lin_space(x2, x3, order);
		std::vector<double> values3 = lin_space(x3, to, order-1);

//		std::vector<double> values = lin_space(from, to, order*2);
		std::vector<double> values;
		for( auto& v : values1 ) values.push_back(v);
		for( auto& v : values2 ) values.push_back(v);
		for( auto& v : values3 ) values.push_back(v);

		std::vector<std::vector<Point>> point_list(N_spins, std::vector<Point>(values.size()));


		for( std::size_t i = 0; i < values.size(); ++i ){
			auto sqrt_s = values[i];
			dummy->mass(sqrt_s);

			auto const& incoming = _diagrams[0]->incoming_particles();
			auto const& outgoing = _diagrams[0]->outgoing_particles();
			Kinematics kin(incoming[0]->mass(), 1, 2);
			auto const q = momentum(incoming[0]->mass(), outgoing[0]->mass(), outgoing[1]->mass());
			kin.incoming(0, four_momentum(0, incoming[0]->mass()));
			kin.outgoing(0, four_momentum(q, outgoing[0]->mass()));
			kin.outgoing(1, four_momentum(-q, outgoing[1]->mass()));

			for( auto& diagram : _diagrams ){
				diagram->reset_spins();
				diagram->reset_indices();
			}

			std::vector<Complex> result(N_spins, 0.);

			for( std::size_t j = 0; j < N_spins; ++j ){
				Complex M{0, 0};
				for( auto& diagram : _diagrams ){
					M += diagram->evaluate_amplitude(kin);
					diagram->iterate_spins();
				}
				point_list[j][i] = Point{sqrt_s, M};
			}
		}

		std::vector<Polynomial> polynomials;
		polynomials.reserve(point_list.size());
		for( auto points : point_list ){
			Polynomial p(order);
			p.fit(points);
			polynomials.push_back(p);
		}
		return polynomials;

	}

	std::vector<Complex> Feynman_Process::decay_M(){
		using namespace Feynumeric::Units;
		validate_diagram_compatibility();

		auto const& incoming = _diagrams[0]->incoming_particles();
		auto const& outgoing = _diagrams[0]->outgoing_particles();

		Kinematics kin(incoming[0]->mass(), 1, 2);

		auto const q = momentum(incoming[0]->mass(), outgoing[0]->mass(), outgoing[1]->mass());

		kin.incoming(0, four_momentum(0, incoming[0]->mass()));
		kin.outgoing(0, four_momentum(q, outgoing[0]->mass()));
		kin.outgoing(1, four_momentum(-q, outgoing[1]->mass()));


		for( auto& diagram : _diagrams ){
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

		std::vector<Complex> result(N_spins, 0.);

		for( std::size_t i = 0; i < N_spins; ++i ){
			Complex M{0, 0};
			for( auto& diagram : _diagrams ){
				result[i] += diagram->evaluate_amplitude(kin);
				diagram->iterate_spins();
			}
		}
		return result;
	}

	double Feynman_Process::decay_1_2(double sqrt_s){
		using namespace Feynumeric::Units;
		validate_diagram_compatibility();

		auto const& incoming = _diagrams[0]->incoming_particles();
		auto const& outgoing = _diagrams[0]->outgoing_particles();

		Kinematics kin(sqrt_s, 1, 2);

		auto const q = momentum(sqrt_s, outgoing[0]->mass(), outgoing[1]->mass());

		kin.incoming(0, four_momentum(0, sqrt_s));
		kin.outgoing(0, four_momentum(q, outgoing[0]->mass()));
		kin.outgoing(1, four_momentum(-q, outgoing[1]->mass()));

		double Ms_squared{0.};

		for( auto& diagram : _diagrams ){
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
//				std::cout << "temp: " << temp << "\n";
				M += temp;
				diagram->iterate_spins();
			}
			Ms_squared += ( M * std::conj(M)).real();
		}
		auto phase_space = 1. / N_polarisations * q / ( 8 * M_PI * sqrt_s * sqrt_s);
//		std::cout << "s: " << sqrt_s << "\n";
//		std::cout << "q: " << q << "\n";
//		std::cout << "N: " << N_polarisations << "\n";
//		std::cout << "phase_space: " << phase_space << "\n";
		return Ms_squared * phase_space;
	}

	double Feynman_Process::partial_decay_1_3(double sqrt_s, double const invariant_mass, double const cos_theta, std::size_t const N_spins, std::size_t const N_polarisations){
		using namespace Feynumeric::Units;

		auto const& incoming = _diagrams[0]->incoming_particles();
		auto const& outgoing = _diagrams[0]->outgoing_particles();

		Kinematics kin(sqrt_s,1, 3);

		auto const q01 = momentum(invariant_mass, outgoing[0]->mass(), outgoing[1]->mass());
		auto const q2 = momentum(sqrt_s, invariant_mass, outgoing[2]->mass());

		//auto const p_out = four_momentum(q2, invariant_mass);

		auto const p_in    = four_momentum(0, sqrt_s);
		auto const p0_rest = four_momentum(q01, outgoing[0]->mass(), cos_theta);
        auto const p1_rest = four_momentum(-q01, outgoing[1]->mass(), cos_theta);
        auto const p01     = four_momentum(q2, invariant_mass);
        auto const p2      = four_momentum(-q2, outgoing[2]->mass());
        auto const p0      = p0_rest.boost(p01);
        auto const p1      = p1_rest.boost(p01);

		kin.incoming(0, p_in);
		kin.outgoing(0, p0);
		kin.outgoing(1, p1);
		kin.outgoing(2, p2);

		// sanity check
		auto p_i0 = kin.incoming(0);
		auto p_f0 = kin.outgoing(0);
		auto p_f1 = kin.outgoing(1);
		auto p_f2 = kin.outgoing(2);
		auto p_ftotal = p_f0 + p_f1 + p_f2;
		auto x1 = p_i0.squared();
		auto x2 = p_ftotal.squared();
		if( !almost_identical(x1, x2, 1.e-3) ){
		    warning(FORMAT("NOT EQUAL {} {}\n", x1, x2));
		}
		if( std::isnan(x1) || std::isnan(x2) ){
		    return 0;
		}

		double Ms_squared{0.};

		for( auto& diagram : _diagrams ){
			diagram->reset_spins();
			diagram->reset_indices();
		}

		for( std::size_t i = 0; i < N_spins; ++i ){
			Complex M{0, 0};
			for( auto& diagram : _diagrams ){
				auto const& temp = diagram->evaluate_amplitude(kin);
				M += temp;
				diagram->iterate_spins();
			}
			Ms_squared += ( M * std::conj(M)).real();
		}

		auto phase_space = 1. / N_polarisations * std::pow(2 * M_PI, -5) * 1. / (16 *  sqrt_s * sqrt_s) * q01 * q2 * 8 * M_PI * M_PI;
		return Ms_squared * phase_space;
	}

	double Feynman_Process::decay_1_3(double sqrt_s){
		using namespace std::placeholders;
		validate_diagram_compatibility();

		auto const M_min = _diagrams[0]->_graph._outgoing[0]->particle()->mass() + _diagrams[0]->_graph._outgoing[1]->particle()->mass();
		auto const M_max = sqrt_s - _diagrams[0]->_graph._outgoing[2]->particle()->mass();
		auto partial = [&](double M){
			auto f = std::bind(&Feynman_Process::partial_decay_1_3, this, sqrt_s, M, _1, _n_spins, _n_polarisations);
			auto r = integrate(f, -0.999, 0.999, 1.e-2);
			return r;
		};
		auto r = integrate(partial, M_min + 1.e-6, M_max-1.e-6, 1.e-2);
		return r;
	}

	double Feynman_Process::decay_width(double sqrt_s){
		if( _diagrams.empty() ){
			critical_error("Process has no diagrams.");
		}
		if( _diagrams[0]->_graph._incoming.size() == 1 ){
			if( _diagrams[0]->_graph._outgoing.size() == 2 ){
				return decay_1_2(sqrt_s);
			}
			if( _diagrams[0]->_graph._outgoing.size() == 3 ){
				return decay_1_3(sqrt_s);
			}
		}
		critical_error("Only single particle decays into a two or three particle final state are implemented.");
	}

	std::map<double, double> Feynman_Process::sigma_table(std::vector<double> const& values, double){
		using namespace Feynumeric::Units;
		using namespace std::placeholders;
        std::map<double, double> result;
        if( _diagrams.empty() ){
            return result;
        }

		std::vector<Feynman_Process> copies;
		copies.reserve(values.size());
		for( std::size_t i = 0; i < values.size(); ++i ){
			copies.emplace_back(*this);
		}

		std::size_t completed{0};
		std::size_t modulo = static_cast<std::size_t>(0.1 * copies.size());
        modulo = modulo == 0? 1 : modulo;

		#pragma omp parallel for
		for( std::size_t i = 0; i < copies.size(); ++i ){
//            std::cout << FORMAT("start {}/{} core: {} {}\n", i, copies.size(), omp_get_thread_num(), values[i]) << std::flush;
			double const sqrt_s = values[i];
			auto f = std::bind(&Feynman_Process::no_check_dsigma_dcos, &copies[i], sqrt_s, _1);
			auto temp = integrate(f, -.999, .999, 0.01);
//            std::cout << "integrate: " << temp << "\n";
			result[sqrt_s] = std::isnan(temp)? 0 : temp;
//			std::cout << FORMAT("end {}/{} core: {} {}\n", i, copies.size(), omp_get_thread_num(), sqrt_s) << std::flush;
            completed++;
            if( completed%modulo == 0 ){
                std::cout << "#" << std::flush;
            }
		}
		std::cout << "\n";
//		std::cout << "End loop\n" << std::flush;
		return result;
	}

	std::map<double, double> Feynman_Process::sigma_table(double start, double end, double delta, double epsilon){
		if( end < start ){
			warning(FORMAT("sigma_table end {} is smaller than start{}.", end, start));
		}
		std::size_t const steps = (end-start) / delta;
		std::vector<double> values(steps + 1);
		for( std::size_t i = 0; i < steps; ++i ){
			values[i] = start + i * delta;
		}
		values[steps] = end;
		return sigma_table(values, epsilon);
	}

	std::map<double, double> Feynman_Process::sigma_table(double start, double end, std::size_t steps, double epsilon){
		if( end < start ){
			warning(FORMAT("sigma_table end {} is smaller than start{}.", end, start));
		}
		std::vector<double> values(steps + 1);
		double const delta = (end-start) / steps;
		for( std::size_t i = 0; i < steps; ++i ){
			values[i] = start + i * delta;
		}
		values[steps] = end;
		return sigma_table(values, epsilon);
	}

	void Feynman_Process::print_sigma_table(std::ostream& out, std::vector<double> const& values, double epsilon){
		auto result = sigma_table(values, epsilon);
        out << "{\n";
        bool first = true;
		for( auto const& [key, value] : result ){
            if( ! first ){
                out << ",";
            }
			out << "{" <<  std::setw(10) << std::fixed << std::setprecision(10) << key << "," << value << "}";
            first = false;
		}
        out << "}\n";
	}

	void Feynman_Process::print_sigma_table(std::ostream& out, double start, double end, double delta, double epsilon){
		auto result = sigma_table(start, end, delta, epsilon);
        out << "{";
        bool first = true;
        for( auto const& [key, value] : result ){
            if( ! first ){
                out << ",";
            }
            out << "{" <<  std::setw(10) << std::fixed << std::setprecision(10) << key << "," << value << "}";
            first = false;
        }
        out << "}";
	}

	void
	Feynman_Process::print_sigma_table(std::ostream& out, double start, double end, std::size_t steps, double epsilon){
		auto result = sigma_table(start, end, steps, epsilon);
        out << "{";
        bool first = true;
        for( auto const& [key, value] : result ){
            if( ! first ){
                out << ",";
            }
            out << "{" <<  std::setw(10) << std::fixed << std::setprecision(10) << key << "," << value << "}";
            first = false;
        }
        out << "}";
	}

    double Feynman_Process::no_check_dsigma_dcos_dM(double sqrt_s, double M, double cos_theta)
    {
        return sqrt_s + M + cos_theta;
    }


	std::vector<std::vector<Polynomial>> Feynman_Process::decay_amplitude1_2(std::vector<double> const& s_values, std::size_t order, bool overwrite_propagator){
		validate_diagram_compatibility();

		std::size_t const N_spins = _diagrams[0]->n_spins();

		std::map<std::size_t, std::vector<std::function<Matrix(Particle::Edge_Ptr, Kinematics const&)>>> virtual_functions;

		for( std::size_t i = 0; i < _diagrams.size() && overwrite_propagator; ++i ){
			for( auto const& vp : _diagrams[i]->_virtual_particles ){
				virtual_functions[i].push_back(vp->feynman_virtual);
				vp->feynman_virtual = [](std::shared_ptr<Graph_Edge> e, Kinematics const& kin){
					return Projector(e, kin);
				};
			}
		}

		auto particle = _diagrams[0]->_incoming_particles[0];
		std::vector<std::vector<Point>> point_list(N_spins, std::vector<Point>(s_values.size()));


		for( std::size_t i = 0; i < s_values.size(); ++i ){
			auto sqrt_s = s_values[i];

			auto const& incoming = _diagrams[0]->incoming_particles();
			auto const& outgoing = _diagrams[0]->outgoing_particles();
			Kinematics kin(sqrt_s, 1, 2);
			auto const q = momentum(sqrt_s, outgoing[0]->mass(), outgoing[1]->mass());
			kin.incoming(0, four_momentum(0, sqrt_s));
			kin.outgoing(0, four_momentum(q, outgoing[0]->mass()));
			kin.outgoing(1, four_momentum(-q, outgoing[1]->mass()));

			for( auto& diagram : _diagrams ){
				diagram->reset_spins();
				diagram->reset_indices();
			}

			std::vector<Complex> result(N_spins, 0.);

			for( std::size_t j = 0; j < N_spins; ++j ){
				Complex M{0, 0};
				for( auto& diagram : _diagrams ){
					M += diagram->evaluate_amplitude(kin);
					diagram->iterate_spins();
				}
				point_list[j][i] = Point{sqrt_s, M};
			}
		}

		std::vector<Polynomial> polynomials;
		polynomials.reserve(point_list.size());
		for( auto points : point_list ){
			Polynomial p(order);
			p.fit(points);
			polynomials.push_back(p);
		}

		for( std::size_t i = 0; i < _diagrams.size() && overwrite_propagator; ++i ){
			for( std::size_t j = 0; j < _diagrams[i]->_virtual_particles.size(); ++j){
				_diagrams[i]->_virtual_particles[j]->feynman_virtual = virtual_functions[i][j];
			}
		}
		return {polynomials};
	}

	std::vector<std::vector<Polynomial>> Feynman_Process::decay_amplitude(std::vector<double> const& s_values, std::size_t order, bool overwrite_propagator){
		if( _diagrams[0]->_outgoing_particles.size() == 2 ){
			return decay_amplitude1_2(s_values, order, overwrite_propagator);
		}
		critical_error("Only 2 outgoing particles are implemented.");
	}

	std::vector<std::vector<Polynomial>> Feynman_Process::scattering_amplitude2_2(std::vector<double> const& s_values, std::vector<std::size_t> const& order, bool overwrite_propagator){
		struct vec2{
			double c, s;
		};
		// sampling
		auto c_values = lin_space(-0.999, 0.999, 2 * order[1] + 1);

		double const fixed_s = *(s_values.begin() + s_values.size()/2);
		double const fixed_c = 0.1337;

		std::array<std::vector<vec2>, 2> sample_points;

		for( auto const& sv : s_values ){
			sample_points[0].push_back({fixed_c, sv});
		}
		for( auto const& cv : c_values ){
			sample_points[1].push_back({cv, fixed_s});
		}


		std::map<std::size_t, std::vector<std::function<Matrix(Particle::Edge_Ptr, Kinematics const&)>>> virtual_functions;

		for( std::size_t i = 0; i < _diagrams.size() && overwrite_propagator; ++i ){
			for( auto const& vp : _diagrams[i]->_virtual_particles ){
				virtual_functions[i].push_back(vp->feynman_virtual);
				vp->feynman_virtual = [](std::shared_ptr<Graph_Edge> e, Kinematics const& kin){
					return Projector(e, kin);
				};
			}
			_diagrams[i]->generate_amplitude();
		}

		std::vector<Feynman_Process> copies(std::max(c_values.size(), s_values.size()), *this);

		std::array<std::vector<std::vector<Point>>, 2> samples{
				std::vector<std::vector<Point>>(_n_spins),
				std::vector<std::vector<Point>>(_n_spins)
		};

		std::vector<std::vector<Polynomial>> result;

		for( std::size_t h = 0; h < sample_points.size(); ++h ){
			#pragma omp parallel for
			for( std::size_t i = 0; i < sample_points[h].size(); ++i ){
				auto const cos_theta = sample_points[h][i].c;
				auto const sqrt_s    = sample_points[h][i].s;

				auto& process = copies[i];
				for( auto& diagram : process._diagrams ){
					diagram->reset_spins();
					diagram->reset_indices();
				}

				Kinematics kin(sqrt_s, 2, 2);

				auto const& incoming = process._diagrams[0]->incoming_particles();
				auto const& outgoing = process._diagrams[0]->outgoing_particles();

				auto qin = momentum(sqrt_s, incoming[0]->mass(), incoming[1]->mass());
				auto qout = momentum(sqrt_s, outgoing[0]->mass(), outgoing[1]->mass());

				kin.incoming(0, four_momentum(qin, incoming[0]->mass(), 1));
				kin.incoming(1, four_momentum(-qin, incoming[1]->mass(), 1));
				kin.outgoing(0, four_momentum(qout, outgoing[0]->mass(), cos_theta));
				kin.outgoing(1, four_momentum(-qout, outgoing[1]->mass(), cos_theta));


				for( std::size_t k = 0; k < process._n_spins; ++k ){
					Complex M{0, 0};
					for( std::size_t j = 0; j < process._diagrams.size(); ++j ){
						auto const& temp = process._diagrams[j]->evaluate_amplitude(kin);
						M += temp;
						process._diagrams[j]->iterate_spins();
					}
					#pragma omp critical
					{
						if( h == 0 ){
							samples[h][k].emplace_back(sqrt_s, M);
						}
						else if( h == 1 ){
							samples[h][k].emplace_back(cos_theta, M);
						}
					}
				}

			}

			std::vector<Polynomial> polynomials;
			for( std::size_t k = 0; k < _n_spins; ++k ){
				Polynomial p(order[h]);
				p.fit(samples[h][k]);
				if( h == 0 ){
					auto rescale = p(fixed_s);
					if( rescale == Complex(0., 0.) ){
						critical_error("rescale factor is zero. Try a different order for sqrt_s polynomial.");
					}
					polynomials.push_back(p/rescale); // rescale
				}else{
					polynomials.push_back(p);
				}
			}
			result.push_back(polynomials);
		}

		for( std::size_t i = 0; i < _diagrams.size() && overwrite_propagator; ++i ){
			for( std::size_t j = 0; j < _diagrams[i]->_virtual_particles.size(); ++j){
				_diagrams[i]->_virtual_particles[j]->feynman_virtual = virtual_functions[i][j];
			}
		}

		return result;
	}

	std::vector<std::vector<Polynomial>> Feynman_Process::scattering_amplitude(std::vector<double> const& s_values, std::vector<std::size_t> const& order, bool overwrite_propagator){
		if( _diagrams[0]->_incoming_particles.size() == 2 && _diagrams[0]->_outgoing_particles.size() == 2){
			return scattering_amplitude2_2(s_values, order, overwrite_propagator);
		}
		critical_error("Only 2->2 scattering processes are implemented.");
	}
}