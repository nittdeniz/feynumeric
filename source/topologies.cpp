#include "topologies.hpp"


namespace Feyncalc
{
    namespace Topology
    {
        [[maybe_unused]] const Graph Double_Wrench(
                vector<int>{0, 1}, vector<int>{2, 3}, vector<int>{4,5},
                vector<array<int, 2>>{{0,2},{1,2},{2,3},{3,4},{3,5}}
                );
    }
}
