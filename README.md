# Feynumeric
A library to numerically evaluate cross sections and decay widths of Feynman diagrams. Highly customizable.

# Requirements / Dependencies

This library uses the Catch2 framework for tests and the fmt library.


# Basics: QED in Feynumeric

To calculate the process `e^+e^- -> µ^+ µ^-` we can use the following script:

```cpp
#include <feynumeric/core.hpp>
#include <feynumeric/qed.hpp> // feynumeric provides the standard QED rules in this header

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
                                    VMP, // A pointer to the Topology_Vertex Manager
	                            {Electron, Positron}, // incoming particles
	                            {Photon}, // virtual particles
	                            {Muon_Minus, Muon_Plus} // outgoing particles
                             );
     
    Feynman_Process pair_production({s_channel}); // add the diagram to the process
	  pair_production.print_dsigma_dcos_table(std::cout, 500._MeV, 0.1); // print to std::cout, sqrt_s is 500 MeV and we go from cos \theta = -1 to 1 in 0.1-sized steps
}

```
# Topologies

Each diagram in Feynumeric is represented by a topology. The topology is agnostic to each edge's particle types. 
We start by drawing a simple graph, labeling the incoming nodes with the letter `i` and an index, starting from `0`. 
Likewise we prepend virtual vertices with a `v` and outgoing vertices with an `o`. 
A standard s-channel scattering process, i.e. annihilation and pair creation would look like this  

````
i0                  o0
 \                 /
  \               /
   \             /
    v0----------v1      
   /             \
  /               \
 /                 \
i1                  o1
````

Translated into code, we add the diagram as pairs in the direction of the particle. I.e. `{"v0", "v1"}` means the 
particle is considered outgoing from the vertex `v0` and incoming at the vertex `v1`. 

````cpp
const Topology s_channel = {{
  {"i0", "v0"},
  {"v0", "v1"},
  {"i1", "v0"},
  {"v1", "o0"},
  {"v1", "o1"}
}};
````

To create your own topology, simply draw the diagram on a piece of paper and label the vertices appropriately.

Interference
===
To calculate interferences, note that all diagrams must have the same particles as `i0` ... `in` and `o0` ... `on`. 


Adding Feynman Rules
====
The class `Feynumeric::Vertex` takes two arguments in its constructor. Firstly, a list of
the particles and their directions w.r.t. the vertex. `Edge_Direction::ANY` is the default argument and
does not necessarily need to be specified. Other options are `Edge_Direction::IN` and `Edge_Direction::OUT`.

````cpp
Feynumeric::Vertex(
    {
        {Electron, Edge_Direction::ANY},
        {Positron, Edge_Direction::ANY},
        {Photon,   Edge_Direction::ANY}
    },
    [](Feynumeric::Kinematics const&, std::vector<std::shared_ptr<Graph_Edge>> const& edges){
        using namespace Feynumeric::Units;
        auto const& photon = edges[2];
        return 1._e * Feynumeric::GAC[*( photon->lorentz_indices()[0] )];
    }
));
````

The second argument is a function that takes a `Feynumeric::Kinematics` object and a vector of `std::shared_ptr<Graph_Edge>`.
`edges` contains the respective edges of the diagram in the same order as the particles defined above. 
In this case, `edge[2]` corresponds to the photon, and we can access its lorentz_indices, momentum or other properties.

Feynman_Process::decay_width
====
To calculate the decay width of a particle, include all contributing diagrams in your feynman process. 
Then call `.decay_width(double mass)` which takes mass of the decaying particle as an argument. 
This way, it is also possible to calculate the decay width for offshell particles.
````cpp
auto diagram_1 = create_diagram(...);
...
auto diagram_n = create_diagram(...);
Feynman_Process decay_process({diagram_1, ..., diagram_n});
auto width = decay_process.decay_width(particle->mass());
````

Feynman_Process::sigma_table
====



Feynman_Process::scattering_amplitude
====
For time consuming operations, _Feynumeric_ offers to export a polynomial fit. 
The fit is done in the following manner (here `c = cos theta` and `s = sqrt_s`).
````
M(c, s) = P(k, s) * Q(c, m) 
````
`P` is a polynomial fit in `s`, evaluated at a constant cosine value `k` and `Q` is a polynomial in `c` on the mass-shell.
The fit points in `s` are weighted 50% within `s +/- width` where a more accurate depiction is necessary.

**Make sure to set the coupling constants to `1` before doing a fit**. This way it is easier to rescale the amplitudes. 