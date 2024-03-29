#include "feynumeric/kinematics.hpp"

namespace Feynumeric
{

	Kinematics::Kinematics(double sqrt_s, std::size_t n_in, std::size_t n_out)
	: _sqrt_s(sqrt_s)
	, _n_in(n_in)
	, _n_out(n_out)
	{
		_momenta.resize(n_in + n_out);
        _angles.resize(n_out-1);
	}

	double Kinematics::sqrt_s() const
	{
		return _sqrt_s;
	}

//    double Kinematics::phi(){
//        double s = _sqrt_s * _sqrt_s;
//        double s2 = s*s;
//        double m1_2 = _momenta[_n_in+0].squared();
//        double m2_2 = _momenta[_n_in+1].squared();
//        auto beta = [](double M1, double M2, double s){
//            return std::sqrt(1- 2*(M1+M2)/s + (M1-M2)*(M1-M2)/(s*s));
//        };
//        if( _n_out == 2 ){
//            return beta(m1_2, m2_2, s)/(32*M_PI * M_PI);
//        }else if( _n_out == 3 ){
//            beta()
//        }
//    }

	Four_Vector const& Kinematics::incoming(std::size_t i) const
	{
		return _momenta[i];
	}

	Four_Vector const& Kinematics::outgoing(std::size_t i) const
	{
		return _momenta[_n_in + i];
	}

	Four_Vector const& Kinematics::momentum(std::size_t i) const
	{
		return _momenta[i];
	}

	void Kinematics::incoming(std::size_t i, Four_Vector const& p)
	{
		_momenta[i] = p;
	}

	void Kinematics::outgoing(std::size_t i, Four_Vector const& p)
	{
		_momenta[_n_in+i] = p;
	}

    void Kinematics::angle(std::size_t i, double a)
    {
        _angles[i] = a;
    }

    double Kinematics::angle(std::size_t i) const
    {
        return _angles[i];
    }


}