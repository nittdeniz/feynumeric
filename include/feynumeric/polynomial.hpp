#ifndef FEYNUMERIC_POLYNOMIAL_HPP
#define FEYNUMERIC_POLYNOMIAL_HPP

#include "complex.hpp"
#include "functions.hpp"
#include "utility.hpp"

#include <array>
#include <vector>

#include <type_traits>

template<typename... T>
struct all_same : std::false_type { };

template<>
struct all_same<> : std::true_type { };

template<typename T>
struct all_same<T> : std::true_type { };

template<typename T, typename... Ts>
struct all_same<T, T, Ts...> : all_same<T, Ts...> { };


namespace Feynumeric
{
	struct Point{
		double x;
		Complex y;

		Point(double x = 0., Complex y = 0.)
		: x(x)
		, y(y)
		{}
		Point(Point const& p)
		: x(p.x)
		, y(p.y){}
		Point& operator=(Point const& p){
			x = p.x;
			y = p.y;
			return *this;
		}
	};

	template <std::size_t N>
	class FPolynomial;

	class Polynomial
	{
		std::size_t n;
		std::size_t order;
		std::vector<Complex> _coefficients;
	public:
		Polynomial(std::size_t order = 0);
		Polynomial(std::vector<Complex> const& coefficients);
		Polynomial(Polynomial const& other);
		void fit(std::vector<Point> const& data);
		Complex integrate(double a, double b);
		std::string to_string(char x) const;

		Polynomial conjugate() const;

		Complex operator()(double x) const;

		Polynomial& operator=(Polynomial const& other);
		Polynomial& operator+=(Polynomial const& other);
		Polynomial& operator-=(Polynomial const& other);
		Polynomial& operator*=(Polynomial const& other);
		Polynomial& operator*=(Complex const& scale);
		Polynomial& operator/=(Complex const& scale);

		void save(std::string const& file_name) const;
		void load(std::string const& file_name);

		friend Polynomial operator+(Polynomial const& lhs, Polynomial const& rhs);
		friend Polynomial operator-(Polynomial const& lhs, Polynomial const& rhs);
		friend Polynomial operator*(Polynomial const& lhs, Polynomial const& rhs);
		friend Polynomial operator*(Polynomial const& lhs, Complex const& rhs);
		friend Polynomial operator*(Complex const& lhs, Polynomial const& rhs);
		friend Polynomial operator/(Polynomial const& lhs, Complex const& rhs);

		template <std::size_t N>
		friend class FPolynomial;
	};


	template <std::size_t N>
	class FPolynomial
	{
		std::size_t n;
		std::size_t order;
		std::vector<func_t<N>> _coefficients;
	public:
		FPolynomial(std::size_t order = 0);
		FPolynomial(std::vector<func_t<1>> const& coefficients);
		FPolynomial(func_t<N> const& f, Polynomial const& p);
		FPolynomial(FPolynomial const& other);
		func_t<N> integrate(double a, double b);
		std::string to_string(char x) const;

		template <typename... T>
		std::enable_if_t<std::conjunction_v<std::is_same<double, T>...>, Complex>
		operator()(double x, T... args) const;

		FPolynomial& operator=(FPolynomial const& other);
		FPolynomial& operator+=(FPolynomial const& other);
		FPolynomial& operator-=(FPolynomial const& other);
		FPolynomial& operator*=(FPolynomial const& other);
		FPolynomial& operator*=(Complex const& scale);
		FPolynomial& operator/=(Complex const& scale);

		template<std::size_t U>
		friend FPolynomial<U> operator+(FPolynomial<U> const& lhs, FPolynomial<U> const& rhs);
		template<std::size_t U>
		friend FPolynomial<U> operator-(FPolynomial<U> const& lhs, FPolynomial<U> const& rhs);
		template<std::size_t U>
		friend FPolynomial<U> operator*(FPolynomial<U> const& lhs, FPolynomial<U> const& rhs);
		template<std::size_t U>
		friend FPolynomial<U> operator*(FPolynomial<U> const& lhs, Complex const& rhs);
		template<std::size_t U>
		friend FPolynomial<U> operator*(Complex const& lhs, FPolynomial<U> const& rhs);
		template<std::size_t U>
		friend FPolynomial<U> operator/(FPolynomial<U> const& lhs, Complex const& rhs);
	};

	template <std::size_t N>
	FPolynomial<N>::FPolynomial(std::vector<func_t<1>> const& coefficients)
			: _coefficients(coefficients)
	{

	}

	template <std::size_t N>
	FPolynomial<N>::FPolynomial(func_t<N> const& f, Polynomial const& p)
	{
		_coefficients.reserve(p.n);
		for( auto const& coef : p._coefficients ){
			auto temp = coef * f;
			_coefficients.push_back(coef * f);
		}
	}

	template <std::size_t N>
	FPolynomial<N>::FPolynomial(FPolynomial<N> const& other)
			: n(other.n)
			  , order(other.order)
			  , _coefficients(other._coefficients)
	{
	}

	template <std::size_t N>
	func_t<N> FPolynomial<N>::integrate(double a, double b){
		func_t<N> result(identity<N>(0.));
		for( std::size_t i = 0; i < _coefficients.size(); ++i ){
			result += _coefficients[i] * identity<N>(int_pow(b, i) - int_pow(a, i));
		}
		return result;
	}

//	template <std::size_t N>
//	std::string FPolynomial<N>::to_string(char x) const{
//		return std::string();
//	}

	template <std::size_t N>
	template <typename... Ts>
	std::enable_if_t<std::conjunction_v<std::is_same<double, Ts>...>, Complex>
	FPolynomial<N>::operator()(double arg, Ts... args) const{
		Complex result{0.};
		for( std::size_t i = 0; i < _coefficients.size(); ++i ){
			result += _coefficients[i](args...) * int_pow(arg, i);
		}
		return result;
	}

	template<std::size_t N>
	FPolynomial<N>::FPolynomial(std::size_t order){

	}

	template<std::size_t N>
	FPolynomial<N>& FPolynomial<N>::operator=(FPolynomial const& other){
		_coefficients = other._coefficients;
		n = other.n;
		order = other.order;
		return *this;
	}
}

#endif