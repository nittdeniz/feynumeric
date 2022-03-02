#include <utility.hpp>
#include "integrate.hpp"

namespace Feynumeric
{
    double integrate(std::function<double(double)> const& f, double const left, double const right, double const epsilon)
    {
	    if( right < left )
	    {
	    	return -integrate(f, right, left);
	    }
    	auto trapezoid = [](double a, double c, double h){return (a+c) * h/2.;};

	    double const delta = right - left;
	    double const half_delta = delta/2.;
    	double const mid  = left + half_delta;
    	double const f_left = f(left);
    	double const f_right = f(right);
    	double const f_mid = f(mid);

    	double const area_total = trapezoid(f_left, f_right, delta);
    	double const area_left  = trapezoid(f_left, f_mid, half_delta);
    	double const area_right = trapezoid(f_mid, f_right, half_delta);

    	if( almost_identical(area_total, area_left + area_right, epsilon) )
	    {
			return area_total;
	    }
    	return integrate(f, left, f_left, mid, f_mid, epsilon) + integrate(f, mid, f_mid, right, f_right, epsilon);
    }

    double integrate(std::function<double(double)> const& f, double left, double f_left, double right, double f_right, double const epsilon)
    {
		if( right < left )
		{
			return -integrate(f, right, f_right, left, f_left, epsilon);
		}
	    auto trapezoid = [](double a, double c, double h){return (a+c) * h/2.;};

	    double const delta = right - left;
	    double const half_delta = delta/2.;
	    double const mid  = left + half_delta;
	    double const f_mid = f(mid);

	    double const area_total = trapezoid(f_left, f_right, delta);
	    double const area_left  = trapezoid(f_left, f_mid, half_delta);
	    double const area_right = trapezoid(f_mid, f_right, half_delta);

	    if( almost_identical(area_total, area_left + area_right, epsilon) )
	    {
		    return area_total;
	    }
	    return integrate(f, left, f_left, mid, f_mid, epsilon) + integrate(f, mid, f_mid, right, f_right, epsilon);
    }
}