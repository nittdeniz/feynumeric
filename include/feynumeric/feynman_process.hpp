#ifndef FEYNUMERIC_FEYNMAN_PROCESS_HPP
#define FEYNUMERIC_FEYNMAN_PROCESS_HPP

#include <initializer_list>
#include <map>
#include <memory>
#include <ostream>
#include <vector>

#include "units.hpp"

namespace Feynumeric
{
	class Feynman_Diagram;
	using Feynman_Diagram_Ptr = std::shared_ptr<Feynman_Diagram>;

	class Feynman_Process
	{
	private:
		std::vector<Feynman_Diagram_Ptr> _diagrams;

		void validate_diagram_compatibility() const;

		double no_check_dsigma_dcos(double sqrt_s, double cos_theta);

		double _conversion_factor = Units::operator""_2barn(1.);

		double decay_1_2();
		double decay_1_3();
		double partial_decay_1_3(double const invariant_mass, double const cos_theta, std::size_t const N_spins, std::size_t const N_polarisations);

		std::size_t _n_spins;
		std::size_t _n_polarisations;

		std::size_t n_spins();
		std::size_t n_polarisations();

	public:
		Feynman_Process(std::initializer_list<Feynman_Diagram_Ptr> list);
		Feynman_Process(std::vector<Feynman_Diagram_Ptr> list);
		Feynman_Process(Feynman_Process const& other);
		void add_diagram(Feynman_Diagram_Ptr diagram);

		void conversion_factor(long double x);

		void print_dsigma_dcos_table(std::ostream& out, double sqrt_s, std::size_t steps);
		void print_dsigma_dcos_table(std::ostream& out, double sqrt_s, double delta);
		void print_dsigma_dcos_table(std::ostream& out, double sqrt_s, std::vector<double>&& values);

		void print_sigma_table(std::ostream& out, std::vector<double> const& values, double epsilon = 1.e-2);
		void print_sigma_table(std::ostream& out, double start, double end, double delta, double epsilon = 1.e-2);
		void print_sigma_table(std::ostream& out, double start, double end, std::size_t steps, double epsilon = 1.e-2);

		std::map<double, double> sigma_table(std::vector<double> const& values, double epsilon = 1.e-2);
		std::map<double, double> sigma_table(double start, double end, double delta, double epsilon = 1.e-2);
		std::map<double, double> sigma_table(double start, double end, std::size_t steps, double epsilon = 1.e-2);



		std::map<double, std::vector<double>> dsigma_dcos_table(double sqrt_s, std::size_t steps);
		std::map<double, std::vector<double>> dsigma_dcos_table(double sqrt_s, double delta);
		std::map<double, std::vector<double>> dsigma_dcos_table(double sqrt_s, std::vector<double>&& values);

		double decay_width();
	};
}
#endif // FEYNUMERIC_FEYNMAN_PROCESS_HPP