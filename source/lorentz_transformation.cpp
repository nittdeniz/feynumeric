#include "feynumeric/lorentz_transformation.hpp"

namespace Feynumeric
{
    Matrix boost(Four_Momentum const& p, Matrix const& a)
    {
        if(
                ( (p.n_cols() == 1 && p.n_rows() == 4) || (p.n_cols() == 4 && p.n_rows() == 1) )
                &&  ( (a.n_cols() == 1 && a.n_rows() == 4) || (a.n_cols() == 4 && a.n_rows() == 1) )
                )
        {
	        // https://en.wikipedia.org/wiki/Lorentz_transformation#Transformation_of_other_quantities
	        Matrix const Z(3, 1, {a.at(1), a.at(2), a.at(3)});
	        Matrix const q(3, 1, {-p.x(), -p.y(), -p.z()}); // formula actually transforms backwards, so we need to change sign
        	Matrix const n(q/std::sqrt(dot3(q, q)));
        	Matrix const beta(q/p.E());
        	double const gamma = 1./std::sqrt(1-dot3(beta, beta).real());
        	auto const& A = a.at(0);
        	auto const Ap = gamma * (A - dot3(beta, Z));
        	auto const Zp = Z + (gamma - 1) * (dot3(Z, n) * n) - gamma * A * beta;

        	return Matrix(a.n_rows(), a.n_cols(), {Ap, Zp.at(0), Zp.at(1), Zp.at(2)});

        }
        critical_error("Invalid boost vector.\n");
    }

	Matrix boost(Four_Momentum const& p, Four_Momentum const& q)
	{
		return boost(p, Matrix(4, 1,{q.E(), q.x(), q.y(), q.z()}));
	}

	Matrix rotateX(Four_Momentum const& p, double cos_theta)
	{
		return Matrix();
	}

	Matrix rotateX(Matrix const& p, double cos_theta)
	{
		return Matrix();
	}

	Matrix rotateY(Four_Momentum const& p, double cos_theta)
	{
		return Matrix();
	}

	Matrix rotateY(Matrix const& p, double cos_theta)
	{
		return Matrix();
	}

	Matrix rotateZ(Four_Momentum const& p, double cos_theta)
	{
		return Matrix();
	}

	Matrix rotateZ(Matrix const& p, double cos_theta)
	{
		return Matrix();
	}

	Matrix project(Four_Momentum const& p, Four_Momentum const& q)
	{

	}
}

