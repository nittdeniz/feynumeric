#ifndef Feynumeric_PROCESS_HPP
#define Feynumeric_PROCESS_HPP

#include <vector>

#include "diagram.hpp"

namespace Feynumeric
{
    using std::vector;

    class Process
    {
    private:
        vector<Diagram> _diagrams;
    public:
        void add_diagrams(vector<Diagram> diagram_list);

        vector<double> dsigma_dcos_theta(double sqrt_s, double cos_theta);
        vector<double> sigma(double sqrt_s);

        double phase_space() const;
    };
}

#endif // Feynumeric_PROCESS_HPP