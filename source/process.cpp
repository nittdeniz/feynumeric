#include "feyncalc/process.hpp"

namespace Feyncalc
{
    void Process::add_diagrams(vector<Diagram<Topology::Double_Wrench>> const&)
    {

    }

    vector<double> Process::dsigma_dcos_theta(double sqrt_s, double cos_theta) const
    {
        return vector<double>({sqrt_s*cos_theta});
    }

    vector<double> Process::sigma(double sqrt_s) const
    {
        return vector<double>({sqrt_s});
    }
}