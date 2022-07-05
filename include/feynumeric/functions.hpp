#pragma once

#include "complex.hpp"

#include <functional>
#include <type_traits>

namespace Feynumeric{
	template <std::size_t, typename T> using alwaysT = T;

	template <typename T, typename Seq> struct func_struct;

	template <typename T, std::size_t... Is>
	struct func_struct<T, std::index_sequence<Is...>> {
		using type = std::function<Complex(alwaysT<Is, T>... floats)>;
	};

	template <std::size_t N>
	struct func_t : func_struct<double, std::make_index_sequence<N>>::type {
		using Base = typename func_struct<double, std::make_index_sequence<N>>::type;
		using Base::Base;
	};

	template <std::size_t N>
	func_t<N> operator+(func_t<N> const& lhs, func_t<N> const& rhs){
		return [lhs, rhs](auto&&... args){
			return lhs(args...) + rhs(args...);
		};
	}

	template <std::size_t N>
	func_t<N> operator-(func_t<N> const& lhs, func_t<N> const& rhs){
		return [lhs, rhs](auto&&... args){
			return lhs(args...) - rhs(args...);
		};
	}

	template <std::size_t N>
	func_t<N> operator*(func_t<N> const& lhs, func_t<N> const& rhs){
		return [lhs, rhs](auto&&... args){
			return lhs(args...) * rhs(args...);
		};
	}

	template <std::size_t N>
	func_t<N> operator*(func_t<N> const& lhs, Complex rhs){
		return [lhs, rhs](auto&&... args){
			return lhs(args...) * rhs;
		};
	}

	template <std::size_t N>
	func_t<N> operator*(Complex lhs, func_t<N> const& rhs){
		return rhs * lhs;
	}

	template <std::size_t N>
	func_t<N> operator/(func_t<N> const& lhs, func_t<N> const& rhs){
		return [lhs, rhs](auto&&... args){
			return lhs(args...) / rhs(args...);
		};
	}

	template <std::size_t N>
	func_t<1> operator>>(func_t<N> const& lhs, std::function<Complex(Complex)> const& rhs){
		return [lhs, rhs](auto&&... args){
			return rhs(lhs(args...));
		};
	}



	inline Complex identity(double x){
		return x;
	}
}