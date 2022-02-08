#ifndef FEYNUMERIC_DIRECTION_HPP
#define FEYNUMERIC_DIRECTION_HPP

#include "messages.hpp"

namespace Feynumeric
{
	enum class Direction : unsigned char
	{
		INCOMING = 0x1,
		OUTGOING = 0x2,
		BOTH     = 0x3,
		VIRTUAL  = 0x4
	};

	constexpr Direction operator|(Direction const& lhs, Direction const& rhs)
	{
		return static_cast<Direction>(static_cast<unsigned char>(lhs) | static_cast<unsigned char>(rhs));
	}

	constexpr Direction operator&(Direction const& lhs, Direction const& rhs)
	{
		return static_cast<Direction>(static_cast<unsigned char>(lhs) & static_cast<unsigned char>(rhs));
	}

	constexpr Direction operator^(Direction const& lhs, Direction const& rhs)
	{
		return static_cast<Direction>(static_cast<unsigned char>(lhs) ^ static_cast<unsigned char>(rhs));
	}

	inline std::string to_string(Direction const& d)
	{
		switch( d )
		{
			case Direction::INCOMING:
				return "incoming";
			case Direction::OUTGOING:
				return "outgoing";
			case Direction::BOTH:
				return "external";
			case Direction::VIRTUAL:
				return "virtual";
			default:
				critical_error("Direction out of bounds. (to_string)");
		}
	}

	inline std::ostream& operator<<(std::ostream& out, Direction const& d)
	{
		return out << to_string(d);
	}
}

#endif // FEYNUMERIC_DIRECTION_HPP