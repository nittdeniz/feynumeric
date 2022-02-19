#ifndef FEYNUMERIC_FEYNMAN_PROCESS_HPP
#define FEYNUMERIC_FEYNMAN_PROCESS_HPP

#include <initializer_list>
#include <map>
#include <memory>
#include <ostream>
#include <vector>

namespace Feynumeric
{
	class Feynman_Diagram;
	using Feynman_Diagram_Ptr = std::shared_ptr<Feynman_Diagram>;

	class Feynman_Process
	{
	private:
		std::vector<Feynman_Diagram_Ptr> _diagrams;

		void validate_diagram_compatibility() const;

	public:
		Feynman_Process(std::initializer_list<Feynman_Diagram_Ptr> list);
		void add_diagram(Feynman_Diagram_Ptr diagram);

		void dsigma_dcos_table(std::ostream& out, double sqrt_s, std::size_t steps);
		void dsigma_dcos_table(std::ostream& out, double sqrt_s, double delta);
		void dsigma_dcos_table(std::ostream& out, double sqrt_s, std::vector<double>&& values);

		std::map<double, std::vector<double>> dsigma_dcos_table(double sqrt_s, std::vector<double>&& values);
	};
}
#endif // FEYNUMERIC_FEYNMAN_PROCESS_HPP