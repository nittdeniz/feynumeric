#include "edge.hpp"
#include "topologies.hpp"
#include <vector>

namespace Feyncalc
{
    namespace Topology
    {
        [[maybe_unused]] const Graph Double_Wrench(
                std::vector<Edge>{
                        {0, 2, Edge::Type::INCOMING},
                        {1, 2, Edge::Type::INCOMING},
                        {2, 3, Edge::Type::VIRTUAL},
                        {3, 4, Edge::Type::OUTGOING},
                        {3, 5, Edge::Type::OUTGOING}
                });

        [[maybe_unused]] const Graph X_Man(
                std::vector<Edge>{
                        {0,2,Edge::Type::INCOMING},
                        {1,3,Edge::Type::INCOMING},
                        {2,3,Edge::Type::VIRTUAL},
                        {2,4,Edge::Type::OUTGOING},
                        {3,5,Edge::Type::OUTGOING}
                }
            );
    }
}
