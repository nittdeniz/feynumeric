#include "feynumeric/momentum.hpp"

namespace Feynumeric
{
    double kallen_lambda(double a, double b, double c)
    {
        return a * a + b * b + c * c - 2 * (a * b + b * c + c * a);
    }

    double momentum(double M, double m1, double m2)
    {
        return std::sqrt(kallen_lambda(M*M, m1*m1, m2*m2)) / (2*M);
    }

    Matrix four_momentum(double mass, double momentum, double cos_theta, double cos_phi)
    {
        const double sin_theta = std::sqrt(1 - cos_theta * cos_theta);
        const double sin_phi   = std::sqrt(1 - cos_phi   * cos_phi);
        return Matrix(4, 1, {std::sqrt(mass*mass + momentum * momentum), momentum * sin_theta * cos_phi, momentum * sin_theta * sin_phi, momentum * cos_theta});
    }

    Three_Momentum::Three_Momentum()
    : Matrix(1, 3, {0,0,0})
    {
    }

    Three_Momentum::Three_Momentum(const Three_Momentum &other)
    : Matrix(1, 3, {other.x(), other.y(), other.z()})
    {
    }

    Three_Momentum &Three_Momentum::operator=(const Three_Momentum &other)
    {
        _data = other._data;
        _rows = other._rows;
        _cols = other._cols;
        return *this;
    }

    double Three_Momentum::x() const
    {
        return _data[0].real();
    }

    double Three_Momentum::y() const
    {
        return _data[1].real();
    }

    double Three_Momentum::z() const
    {
        return _data[2].real();
    }

    void Three_Momentum::x(double x)
    {
        _data[0] = x;
    }

    void Three_Momentum::y(double y)
    {
        _data[1] = y;
    }

    void Three_Momentum::z(double z)
    {
        _data[2] = z;
    }

    double Three_Momentum::dot() const
    {
        return (_data[0] * _data[0] + _data[1] * _data[1] + _data[2] * _data[2]).real();
    }

    Four_Momentum::Four_Momentum()
    : Matrix(1, 4, {0,0,0,0})
    {

    }

	Four_Momentum::Four_Momentum(double q, double m, double cos_theta, double cos_phi)
	: Matrix(1, 4, {0,0,0,0})
	{
    	const double sin_theta = std::sqrt(1 - cos_theta * cos_theta);
    	const double sin_phi   = std::sqrt(1 - cos_phi * cos_phi);
		_data[0] = std::sqrt(q*q + m*m);
		_data[1] = q * cos_phi * sin_theta;
		_data[2] = q * sin_phi * sin_theta;
		_data[3] = q * cos_theta;
	}

	Four_Momentum::Four_Momentum(const Four_Momentum &other)
    : Matrix(1, 4, {other.E(), other.x(), other.y(), other.z()})
    {

    }

    Four_Momentum &Four_Momentum::operator=(const Four_Momentum &other)
    {
        _data = other._data;
        _rows = other._rows;
        _cols = other._cols;
        return *this;
    }

    double Four_Momentum::E() const
    {
        return _data[0].real();
    }

    double Four_Momentum::x() const
    {
        return _data[1].real();
    }

    double Four_Momentum::y() const
    {
        return _data[2].real();
    }

    double Four_Momentum::z() const
    {
        return _data[3].real();
    }

    void Four_Momentum::E(double E)
    {
        _data[0] = E;
    }

    void Four_Momentum::x(double x)
    {
        _data[1] = x;
    }

    void Four_Momentum::y(double y)
    {
        _data[2] = y;
    }

    void Four_Momentum::z(double z)
    {
        _data[3] = z;
    }

    double Four_Momentum::squared() const
    {
        return (
                  _data[0] * std::conj(_data[0])
                - _data[1] * std::conj(_data[1])
                - _data[2] * std::conj(_data[2])
                - _data[3] * std::conj(_data[3])
        ).real();
    }

	Four_Momentum operator+(Four_Momentum const& lhs, Four_Momentum const& rhs)
	{
		Four_Momentum result(lhs);
		result._data[0] += rhs._data[0];
		result._data[1] += rhs._data[1];
		result._data[2] += rhs._data[2];
		result._data[3] += rhs._data[3];
		return result;
	}

	Four_Momentum operator-(Four_Momentum const& lhs, Four_Momentum const& rhs)
	{
		Four_Momentum result(lhs);
		result._data[0] -= rhs._data[0];
		result._data[1] -= rhs._data[1];
		result._data[2] -= rhs._data[2];
		result._data[3] -= rhs._data[3];
		return result;
	}
}