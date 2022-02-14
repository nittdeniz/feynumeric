#ifndef FEYNUMERIC_FEYNMAN_PROCESS_HPP
#define FEYNUMERIC_FEYNMAN_PROCESS_HPP

#include <initializer_list>
#include <ostream>
#include <vector>

namespace Feynumeric
{
	class Feynman_Diagram;

	class Feynman_Process
	{
	private:
		std::vector<Feynman_Diagram*> _diagrams;

		void validate_diagram_compatibility() const;

	public:
		Feynman_Process(std::initializer_list<Feynman_Diagram*> list);
		void add_diagram(Feynman_Diagram* diagram);

		void dsigma_dcos_table(std::ostream& out, double sqrt_s, std::size_t steps);
		void dsigma_dcos_table(std::ostream& out, double sqrt_s, double delta);
		void dsigma_dcos_table(std::ostream& out, double sqrt_s, std::vector<double>&& values);
	};
}
#endif // FEYNUMERIC_FEYNMAN_PROCESS_HPP