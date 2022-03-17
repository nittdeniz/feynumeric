#include "particle.hpp"

namespace Feynumeric
{
    Particle::Particle(std::string&& name, Type type, double mass, double width, double charge, double spin)
    : _name(std::move(name))
    , _type(type)
    , _mass(mass)
    , _charge(charge)
    , _spin(spin, spin, mass==0)
    , _width(width)
    , _baryon_number(0)
    , _lepton_number(0)
    {
    	int sign;
    	if( is_set(type, Type::TrueParticle) ) sign = 1;
    	else if( is_set(type, Type::AntiParticle) ) sign = -1;
    	else if( is_set(type, Type::NeutralParticle) ) sign = 0;
    	else{ warning("Particle type does not contain valid type (Neither True, Anti nor Neutral)"); sign = -666;}

    	if( is_set(type, Type::Baryon) )
	    {
			_baryon_number = sign * 1.;
	    }
    	if( is_set(type, Type::Lepton) )
	    {
    		_lepton_number = sign * 1.;
	    }

    	if( is_set(type, Type::Fermion) && !_spin.is_half_odd_integer() )
	    {
    		critical_error(FORMAT("Particle {} is classified as a fermion but has spin {}.", _name, _spin.j()));
	    }
    	if( type & Type::Boson && _spin.is_half_odd_integer() )
	    {
    		critical_error(FORMAT("Particle {} is classified as a boson but has spin {}.", _name, _spin.j()));
	    }
    }

	Particle::Particle(Particle const& copy)
	: _name(copy._name)
	, _type(copy._type)
	, _mass(copy._mass)
	, _charge(copy._charge)
	, _spin(copy._spin)
	, _isospin(copy._isospin)
	, _parity(copy._parity)
	, _width_function(copy._width_function)
	, _width(copy._width)
	, _baryon_number(copy._baryon_number)
	, _lepton_number(copy._lepton_number)
	, _user_data(copy._user_data)
	{

	}

	Particle& Particle::operator=(Particle const& copy)
	{
		_name           = copy._name;
		_type           = copy._type;
		_mass           = copy._mass;
		_charge         = copy._charge;
		_spin           = copy._spin;
		_isospin        = copy._isospin;
		_parity         = copy._parity;
		_width_function = copy._width_function;
		_width          = copy._width;
		_baryon_number  = copy._baryon_number;
		_lepton_number  = copy._lepton_number;
		_user_data      = copy._user_data;
		return *this;
	}

	double Particle::width(double p2) const
	{
    	if( !_width_function )
	    {
    		return _width;
	    }
		return _width_function(p2);
	}

	void Particle::isospin(Angular_Momentum const& i)
	{
		_isospin = i;
	}

	double Particle::width() const
	{
		return _width;
	}

	void Particle::width(std::function<double(double)> f)
	{
		_width_function = f;
	}

	std::any Particle::user_data(std::string key) const
    {
    	try{
		    return _user_data.at(key);
	    }
    	catch( std::out_of_range const& e)
	    {
    		critical_error(FORMAT("Key <{}> is not set for particle {}\n", key, _name));
	    }
    }

	void Particle::parity(int p)
	{
		_parity = p;
	}

	int Particle::parity() const
	{
		return _parity;
	}

	std::string Particle::name() const
    {
        return _name;
    }

    double Particle::mass() const
    {
        return _mass;
    }

	double Particle::baryon_number() const
	{
		return _baryon_number;
	}

	double Particle::lepton_number() const
	{
		return _lepton_number;
	}

	void Particle::baryon_number(double n)
	{
		_baryon_number = n;
	}

	void Particle::lepton_number(double n)
	{
		_lepton_number = n;
	}

	double Particle::charge() const
    {
        return _charge;
    }

    Angular_Momentum Particle::spin() const
    {
        return _spin;
    }

    Angular_Momentum Particle::isospin() const
    {
        return _isospin;
    }

    void Particle::user_data(std::string key, std::any data)
    {
        _user_data[key] = data;
    }

    bool Particle::is_fermion() const
    {
    	return (_type & Type::Fermion) == static_cast<uint64_t>(Type::Fermion);
    }

    bool Particle::is_true_fermion() const
    {
        return (_type & Type::TrueFermion) == static_cast<uint64_t>(Type::TrueFermion);
    }

    bool Particle::is_anti_fermion() const
    {
        return (_type & Type::AntiFermion) == static_cast<uint64_t>(Type::AntiFermion);
    }

    unsigned int Particle::n_lorentz_indices() const
    {
        return static_cast<unsigned int>(_spin.j());
    }

    std::ostream &operator<<(std::ostream& out, const Particle &p)
    {
        out << p._name;
        return out;
    }

    bool is_fermion(const Particle &particle)
    {
        return particle.is_true_fermion();
    }

    bool is_anti_fermion(const Particle &particle)
    {
        return particle.is_anti_fermion();
    }
}
