#ifndef Feynumeric_PROCESS_HPP
#define Feynumeric_PROCESS_HPP

#include <vector>

#include "feynman_diagram.hpp"

namespace Feynumeric
{
    class Process
    {
    private:
        std::vector<Feynman_Diagram> _diagrams;
    public:
        void add_diagrams(std::vector<Feynman_Diagram> diagram_list);

        std::vector<double> dsigma_dcos_theta(double sqrt_s, double cos_theta);
        std::vector<double> sigma(double sqrt_s);

        double phase_space() const;
    };
}

#endif // Feynumeric_PROCESS_HPP