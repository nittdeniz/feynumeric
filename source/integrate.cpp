#include <random>
#include <utility.hpp>
#include "integrate.hpp"

namespace Feynumeric
{
	double
	integrate(std::function<double(double)> const& f, double const left, double const right, double const epsilon)
    {
        if( right < left )
        {
            return -integrate(f, right, left);
        }
        auto simpson = [](double f_a, double f_mid, double f_b, double d)
        { return d / 6. * (f_a + 4 * f_mid + f_b); };

        double const a = left;
        double const b = a + (right - left) / 1.99981274; // do not go to exactly the middle due to symmetry reasons
        double const c = right;
        double const ab = a + (b - a) / 2.;
        double const bc = b + (c - b) / 2.;

        double const f_a = f(a);
        double const f_ab = f(ab);
        double const f_b = f(b);
        double const f_bc = f(bc);
        double const f_c = f(c);

        double const area_left = simpson(f_a, f_ab, f_b, b - a);
        double const area_right = simpson(f_b, f_bc, f_c, c - b);
        double const area_total = simpson(f_a, f_b, f_c, c - a);

        if( almost_identical(area_total, area_left + area_right, epsilon) || std::isnan(area_total))
        {
            return area_total;
        }
        return integrate(f, a, f_a, ab, f_ab, b, f_b, epsilon) + integrate(f, b, f_b, bc, f_bc, c, f_c, epsilon);
    }

    double integrate(const std::function<double(double)> &f, const double a, const double f_a, const double b, const double f_b,
                                 const double c, const double f_c, const double epsilon)
    {
        auto simpson = [](double f_a, double f_mid, double f_b, double d)
        { return d / 6. * (f_a + 4 * f_mid + f_b); };

        double const ab = a + (b - a) / 2.;
        double const bc = b + (c - b) / 2.;

        double const f_ab = f(ab);
        double const f_bc = f(bc);

        double const area_left = simpson(f_a, f_ab, f_b, b - a);
        double const area_right = simpson(f_b, f_bc, f_c, c - b);
        double const area_total = simpson(f_a, f_b, f_c, c - a);

        if( almost_identical(area_total, area_left + area_right, epsilon) || std::isnan(area_total))
        {
            return area_total;
        }
        return integrate(f, a, f_a, ab, f_ab, b, f_b, epsilon) + integrate(f, b, f_b, bc, f_bc, c, f_c, epsilon);
    }


//
//
//		auto trapezoid = [](double a, double c, double h){ return ( a + c ) * h / 2.; };
//
//		double const delta = right - left;
//		double const half_delta = delta / 2.;
//		double const mid = left + half_delta;
//		double const f_left = f(left);
//		double const f_right = f(right);
//		double const f_mid = f(mid);
//
//		double const area_total = trapezoid(f_left, f_right, delta);
//		double const area_left = trapezoid(f_left, f_mid, half_delta);
//		double const area_right = trapezoid(f_mid, f_right, half_delta);

//		if( almost_identical(area_total, area_left + area_right, epsilon) || std::isnan(area_total) ){
//			return area_total;
//		}
//		return integrate(f, left, f_left, mid, f_mid, epsilon) + integrate(f, mid, f_mid, right, f_right, epsilon);
//	}

	double integrate(std::function<double(double)> const& f, double left, double f_left, double right, double f_right,
	                 double const epsilon){
		if( right < left ){
			return -integrate(f, right, f_right, left, f_left, epsilon);
		}
		auto trapezoid = [](double a, double c, double h){ return ( a + c ) * h / 2.; };

		double const delta = right - left;
		double const half_delta = delta / 2.;
		double const mid = left + half_delta;
		double const f_mid = f(mid);

		double const area_total = trapezoid(f_left, f_right, delta);
		double const area_left = trapezoid(f_left, f_mid, half_delta);
		double const area_right = trapezoid(f_mid, f_right, half_delta);

		if( almost_identical(area_total, area_left + area_right, epsilon) || std::isnan(area_total)){
			return area_total;
		}
		return integrate(f, left, f_left, mid, f_mid, epsilon) + integrate(f, mid, f_mid, right, f_right, epsilon);
	}
}