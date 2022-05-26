#include "particle.hpp"
#include "particle_manager.hpp"
#include "units.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>

namespace Feynumeric
{
	void Particle_Manager::parse_file(std::string const& file_name)
	{
		auto parse_fraction = [](std::string&& str)
		{
			auto pos = std::find(str.begin(), str.end(), '/');
			if( pos == str.end() )
			{
				return std::stod(str);

			}
			auto numerator   = std::stod(std::string(str.begin(), pos));
			auto denominator = std::stod(std::string(pos+1, str.end()));
			return numerator / denominator;
		};

		auto parse_unit = [](std::string&& str)
		{
			using namespace Units;
			std::size_t pos = 0;
			for( ; pos < str.size(); ++pos )
			{
				char const& c = str[pos];
				if( ('A' <= c && c <= 'Z') || ('a' <= c && c <= 'z')  || c == '%')
				{
					break;
				}
			}
			double value = std::stod(std::string(str.begin(), str.begin()+pos));
			std::string unit  = std::string(str.begin()+pos, str.end());
			std::transform(unit.begin(), unit.end(), unit.begin(), [](unsigned char c){ return std::tolower(c);});
			if( unit == "ev"  ) return static_cast<double>(value * 1._eV);
			if( unit == "kev" ) return static_cast<double>(value * 1._keV);
			if( unit == "mev" ) return static_cast<double>(value * 1._MeV);
			if( unit == "gev" ) return static_cast<double>(value * 1._GeV);
			if( unit == "tev" ) return static_cast<double>(value * 1._TeV);
			if( unit == "%" ) return static_cast<double>(value * 1._percent);
			if( unit == "" )  return static_cast<double>(value);
			critical_error(FORMAT("Could not read unit. Supported: eV, keV, MeV, GeV, TeV, %. Found: {}", unit));
		};

		enum class STATE{
			READ_DECLARATION,
			READ_DECLARATION_NAME,
			READ_PARAMETERS,
			END,
			ERROR
		};

		std::ifstream in(file_name);
		if( !in.is_open() )
		{
			critical_error(FORMAT("Could not open file: {}", file_name));
		}
		char c;
		STATE state = STATE::READ_DECLARATION;
		Particle particle;
		bool quit = false;
		while( !quit ){
			std::string input;
			quit = !static_cast<bool>((in >> std::noskipws >> c));
			if( quit ) state = STATE::END;
			else{
				if( !std::isspace(c))
					in.putback(c);
				else{
					continue;
				}
				if( state == STATE::READ_PARAMETERS && c == '@' ){
					state = STATE::END;
				} else{
					in >> std::skipws >> input;
					if( input.empty() ){
						continue;
					}
				}
			}
            if( c == '#' ){ // ignore the comment
                std::getline(in, input);
                continue;
            }
			switch( state ){
				case STATE::READ_DECLARATION:{
					if( input[0] != '@' ){
						critical_error(FORMAT("Parse error: Expected declaration sign (@). Found: {}", c));
					}
					particle = Particle(); // reset
					std::string_view declaration(input.begin() + 1, input.end());
					if( declaration == "group" ){
						particle.is_group(true);
					} else if( declaration == "particle" ){
						particle.is_group(false);
					} else{
						critical_error(FORMAT("Parse error: Declaration must be group or particle. Found: {}", declaration));
					}
					state = STATE::READ_DECLARATION_NAME;
					break;
				}
				case STATE::READ_DECLARATION_NAME:{
					particle._name = input;
					state = STATE::READ_PARAMETERS;
					break;
				}
				case STATE::READ_PARAMETERS:{
					std::string temp;
					do{
						in >> std::skipws >> temp;
					}while( temp.empty() );
					if( input == "mass" ){
						particle._mass = parse_unit(std::move(temp));
					} else if( input == "type" ){
						static const std::map<std::string, Particle::Type> types{{
								{"TrueParticle", Particle::Type::TrueParticle},
								{"AntiParticle", Particle::Type::AntiParticle},
								{"NeutralParticle", Particle::Type::NeutralParticle},
								{"Boson", Particle::Type::Boson},
								{"TrueBoson", Particle::Type::TrueBoson},
								{"AntiBoson", Particle::Type::AntiBoson},
								{"NeutralBoson", Particle::Type::NeutralBoson},
								{"Fermion", Particle::Type::Fermion},
								{"TrueFermion", Particle::Type::TrueFermion},
								{"AntiFermion", Particle::Type::AntiFermion},
								{"NeutralFermion", Particle::Type::NeutralFermion},
								{"Lepton", Particle::Type::Lepton},
								{"TrueLepton", Particle::Type::TrueLepton},
								{"AntiLepton", Particle::Type::AntiLepton},
								{"Meson", Particle::Type::Meson},
								{"TrueMeson", Particle::Type::TrueMeson},
								{"AntiMeson", Particle::Type::AntiMeson},
								{"NeutralMeson", Particle::Type::NeutralMeson},
								{"Baryon", Particle::Type::Baryon},
								{"TrueBaryon", Particle::Type::TrueBaryon},
								{"AntiBaryon", Particle::Type::AntiBaryon},
						}};
						if( types.contains(temp) ){
							particle.set_type(types.at(temp));
						}else{
							critical_error(FORMAT("Unknown type: {}.", temp));
						}
					} else if( input == "spin" ){
						double j = parse_fraction(std::move(temp));
						particle._spin = Angular_Momentum(j, j);
					} else if( input == "isospin" ){
						double j = parse_fraction(std::move(temp));
						in >> temp;
						double m = parse_fraction(std::move(temp));
						particle._isospin = Angular_Momentum(j, m);
					} else if( input == "width" ){
						particle._width = parse_unit(std::move(temp));
					} else if( input == "parity" ){
						particle._parity = std::stoi(temp);
					}else if( input == "charge" ){
						particle._charge = std::stod(temp);
					} else if( input == "group" ){
						if( !_particles.contains(temp) ){
							critical_error(FORMAT("Group {} must be defined before using it.", temp));
						}
						particle.copy_parameters(*_particles[temp]);
						particle._parent = _particles[temp];
					} else if( input == "baryonnumber" ){
						particle._baryon_number = std::stod(temp);
					} else if( input == "leptonnumber" ){
						particle._lepton_number = std::stod(temp);
					}
					else{
						// user data
						particle.user_data(input, parse_unit(std::move(temp)));
						//warning(FORMAT("Parameter {} ignored.", input));
					}
					break;
				}
				case STATE::END:{
					if( _particles.contains(particle._name) ){
						warning(FORMAT("Overwriting particle {}.", particle._name));
					}
					_particles[particle._name] = std::make_shared<Particle>(particle);
					state = STATE::READ_DECLARATION;
					break;
				}
				case STATE::ERROR:
				default:{
					break;
				}
			}
		}
	}

	Particle_Ptr const& Particle_Manager::get(std::string const& key) const{
		if( _particles.contains(key) ){
			return _particles.at(key);
		}
		critical_error(FORMAT("Particle {} does not exist in Particle Manager.", key));
	}

	bool Particle_Manager::exists(std::string const& key) const{
		return _particles.contains(key);
	}

	Particle_Manager::Particle_Manager(std::string const& file_name)
	{
		parse_file(file_name);
	}

	Particle_Manager::Particle_Manager(Particle_Manager const&)
	{

	}

	Particle_Manager& Particle_Manager::operator=(Particle_Manager const&)
	{
		return *this;
	}

	Particle_Ptr const& Particle_Manager::operator[](std::string const& key) const
	{
		try{
			return _particles.at(key);
		} catch( std::out_of_range const& e ){
			critical_error(FORMAT("Particle {} is not defined.", key));
		}
	}

	Particle_Ptr& Particle_Manager::operator[](std::string const& key)
	{
		return _particles[key];
	}

	std::map<std::string, Particle_Ptr>::iterator Particle_Manager::begin(){
		return _particles.begin();
	}

	std::map<std::string, Particle_Ptr>::iterator Particle_Manager::end(){
		return _particles.end();
	}

	std::map<std::string, Particle_Ptr>::const_iterator Particle_Manager::cbegin() const{
		return _particles.cbegin();
	}

	std::map<std::string, Particle_Ptr>::const_iterator Particle_Manager::cend() const{
		return _particles.cend();
	}
}