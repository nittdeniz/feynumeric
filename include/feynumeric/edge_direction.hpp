#ifndef FEYNUMERIC_EDGE_DIRECTION_HPP
#define FEYNUMERIC_EDGE_DIRECTION_HPP

namespace Feynumeric
{
	enum class Edge_Direction : unsigned char
	{
		IN = 0x01,
		OUT = 0x02,
		ANY = 0x03
	};
}

#endif