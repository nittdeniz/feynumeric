#include "feynumeric/process.hpp"

namespace Feynumeric
{
    void Process::add_diagrams(vector<Diagram> diagrams)
    {
        _diagrams.reserve(_diagrams.size() + diagrams.size());
        for( auto& diagram : diagrams )
        {
            diagram.generate_amplitude();
            _diagrams.emplace_back(diagram);
        }

    }

    /**
     * @brief
     * @param sqrt_s
     * @param cos_theta
     * @return an array with components {coherent_sum, diagram_1, diagram_2, ...}
     */
    vector<double> Process::dsigma_dcos_theta(double sqrt_s, double cos_theta)
    {
        const double PS = phase_space();

        vector<double> result(_diagrams.size()+1, 0);
        Complex coherent_sum;
        int i = 1;
        for( auto& diagram : _diagrams )
        {
            const auto M = diagram.calculate_amplitude(sqrt_s, cos_theta);
            result[i++] = PS * (M * std::conj(M)).real();
            coherent_sum += M;
        }
        result[0] = PS * (coherent_sum * std::conj(coherent_sum)).real();
        return result;
    }

    vector<double> Process::sigma(double sqrt_s)
    {
        return vector<double>({sqrt_s});
    }

    double Process::phase_space() const
    {
        return 1;
    }
}