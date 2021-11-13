#ifndef FEYNCALC_PROCESS_HPP
#define FEYNCALC_PROCESS_HPP

#include <vector>

#include "diagram.hpp"

namespace Feyncalc
{
    using std::vector;

    class Process
    {
    private:
        vector<Diagram> _diagrams;
    public:
        void add_diagrams(vector<Diagram> diagram_list);

        vector<double> dsigma_dcos_theta(double sqrt_s, double cos_theta) const;
        vector<double> sigma(double sqrt_s) const;

        double phase_space() const;
    };
}

#endif // FEYNCALC_PROCESS_HPP