#include "direction.hpp"
#include "topologies.hpp"

namespace Feynumeric
{
	const Topology Decay_1_to_2({
			                     {0, 1, Direction::INCOMING},
			                     {1, 2, Direction::OUTGOING},
			                     {1, 3, Direction::OUTGOING}
	});

	const Topology Decay_1_to_3({
			                            {0, 1, Direction::INCOMING},
			                            {1, 2, Direction::VIRTUAL},
			                            {2, 3, Direction::OUTGOING},
			                            {2, 4, Direction::OUTGOING},
			                            {1, 5, Direction::OUTGOING}
	                            });

	const Topology Scattering_Vertical_2_to_2({
			                     {0, 2, Direction::INCOMING},
			                     {1, 3, Direction::INCOMING},
			                     {2, 3, Direction::VIRTUAL},
			                     {2, 4, Direction::OUTGOING},
			                     {3,5, Direction::OUTGOING}
	                     });

	const Topology Scattering_Horizontal_2_to_2({
			                             {0, 2, Direction::INCOMING},
			                             {1,2, Direction::INCOMING},
			                             {2,3, Direction::VIRTUAL},
			                             {3, 4, Direction::OUTGOING},
			                             {3, 5, Direction::OUTGOING}
	                             });
}
