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
    public:
        void add_diagrams(vector<Diagram<Topology::Double_Wrench>> const& diagram_list);

        vector<double> dsigma_dcos_theta(double sqrt_s, double cos_theta) const;
        vector<double> sigma(double sqrt_s) const;
    };
}

#endif // FEYNCALC_PROCESS_HPP