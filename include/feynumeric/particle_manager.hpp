#ifndef FEYNUMERIC_PARTICLE_MANAGER_HPP
#define FEYNUMERIC_PARTICLE_MANAGER_HPP

#include <map>
#include <string>

#include "particle.hpp"

namespace Feynumeric
{
	class Particle_Manager{
	private:
		void parse_file(std::string const& file_name);
		std::map<std::string, Particle_Ptr> _particles;
	public:
		Particle_Manager(std::string const& file_name);
		Particle_Manager(Particle_Manager const& copy);
		Particle_Manager& operator=(Particle_Manager const& copy);

		Particle_Ptr get(std::string&& key) const;

		Particle_Ptr const& operator[](std::string const& key) const;
		Particle_Ptr& operator[](std::string const& key);
	};
}

#endif