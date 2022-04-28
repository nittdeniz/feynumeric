# Feynumeric
A library to numerically evaluate cross sections and decay widths of Feynman diagrams. Highly customizable.

# Requirements / Dependencies

This library uses the Catch2 framework for tests and the fmt library.

# Pair Production | A First Example

To calculate the process `e^+e^- -> µ^+ µ^-` we can use the following script:

```cpp
#include <feynumeric/core.hpp>
#include <feynumeric/qed.hpp> // feyncalc provides the standard QED rules in this header

int main()
{
    using namespace Feynumeric;
    using namespace Feynumeric::Units;
    using namespace QED;
    
    // initialize global variables
    init_particles();
    init_vertices();
    
    auto s_channel = create_diagram("s_channel", // Name of the diagram (can be anything)
                                    s_channel, // Topology (see more below)
                                    VMP, // A pointer to the Vertex Manager
	                            {Electron, Positron}, // incoming particles
	                            {Photon}, // virtual particles
	                            {Muon_Minus, Muon_Plus} // outgoing particles
                             );
     
    Feynman_Process pair_production({s_channel}); // add the diagram to the process
	  pair_production.print_dsigma_dcos_table(std::cout, 500._MeV, 0.1); // print to std::cout, sqrt_s is 500 MeV and we go from cos \theta = -1 to 1 in 0.1-sized steps
}

```
