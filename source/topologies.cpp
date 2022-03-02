#include "direction.hpp"
#include "topologies.hpp"

namespace Feynumeric
{
	const Topology Wrench({
			                     {0, 1, Direction::INCOMING},
			                     {1, 2, Direction::OUTGOING},
			                     {1, 3, Direction::OUTGOING}
	});

	const Topology X_Man({
			                     {0, 2, Direction::INCOMING},
			                     {1, 3, Direction::INCOMING},
			                     {2, 3, Direction::VIRTUAL},
			                     {2, 4, Direction::OUTGOING},
			                     {3,5, Direction::OUTGOING}
	                     });

	const Topology Double_Wrench({
			                             {0, 2, Direction::INCOMING},
			                             {1,2, Direction::INCOMING},
			                             {2,3, Direction::VIRTUAL},
			                             {3, 4, Direction::OUTGOING},
			                             {3, 5, Direction::OUTGOING}
	                             });

	const Topology Scissors( {
			                         {0, 2, Direction::INCOMING},
			                         {1, 3, Direction::INCOMING},
			                         {2, 3, Direction::VIRTUAL},
			                         {3, 4, Direction::OUTGOING},
			                         {2, 5, Direction::OUTGOING}
	                         });
}
