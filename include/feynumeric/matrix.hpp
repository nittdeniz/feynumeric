#ifndef Feynumeric_MATRIX_HPP
#define Feynumeric_MATRIX_HPP

#include "complex.hpp"
#include <functional>
#include <vector>

namespace Feynumeric
{
    class Matrix
    {
    protected:
        std::size_t _rows, _cols;
        std::vector<Complex> _data;
        std::size_t index(size_t row, std::size_t col) const;
    public:
        class dimension_exception : public std::exception
        {
        public:
            dimension_exception(){}
        };

        Matrix(size_t rows = 0, std::size_t cols = 0, std::vector<Complex>&& data = {});
        Matrix(size_t rows, std::size_t cols, Complex diagonal);
        Matrix(Matrix const& other);
        Matrix(Matrix&& other);
        Matrix& operator=(Matrix const& other);
        Matrix& operator=(Matrix&& other);

        std::size_t n_rows() const;
        std::size_t n_cols() const;
        std::size_t n_total() const;

        Complex const& at(size_t i) const;
        Complex& operator[](size_t i);
        Complex& operator()(size_t i, std::size_t j);
        Complex at(size_t i, std::size_t j) const;

        Matrix& operator+=(Matrix const& rhs);
        Matrix& operator-=(Matrix const& rhs);
        Matrix& operator*=(Matrix const& rhs);
        Matrix& operator*=(Complex const& rhs);
        Matrix& operator/=(Complex const& rhs);

        Matrix& apply(std::function<Complex(Complex const&)> const& f);

        Matrix T() const;

        Complex try_as_complex() const;

        template<typename T>
        friend Matrix operator+(Matrix const& lhs, T const& rhs);
        template<typename T>
        friend Matrix operator+(T const& lhs, Matrix const& rhs);
        friend Matrix operator+(Matrix const& lhs, Matrix const& rhs);

        friend Matrix operator-(Matrix const& lhs);

	    template<typename T>
	    friend Matrix operator-(Matrix const& lhs, T const& rhs);
	    template<typename T>
	    friend Matrix operator-(T const& lhs, Matrix const& rhs);
	    friend Matrix operator-(Matrix const& lhs, Matrix const& rhs);

	    template<typename T>
	    friend Matrix operator*(Matrix const& lhs, T const& rhs);
	    template<typename T>
	    friend Matrix operator*(T const& lhs, Matrix const& rhs);
	    friend Matrix operator*(Matrix const& lhs, Matrix const& rhs);

	    template<typename T>
	    friend Matrix operator/(Matrix const& lhs, T const& rhs);

        friend bool operator==(Matrix const& lhs, Matrix const& rhs);
        friend bool operator!=(Matrix const& lhs, Matrix const& rhs);

        friend bool same_dimension(Matrix const& lhs, Matrix const& rhs);
        friend std::ostream& operator<<(std::ostream& out, Matrix const& matrix);
    };

	template<typename T>
	Matrix operator+(Matrix const& lhs, T const& rhs){
		return lhs + Matrix(lhs._cols, lhs._rows, static_cast<Complex>(rhs));
	}
	template<typename T>
	Matrix operator+(T const& lhs, Matrix const& rhs){
		return rhs + lhs;
	}
	Matrix operator+(Matrix const& lhs, Matrix const& rhs);

	Matrix operator-(Matrix const& lhs);

	template<typename T>
	Matrix operator-(Matrix const& lhs, T const& rhs)
	{
		return lhs - Matrix(lhs._cols, lhs._rows, static_cast<Complex>(rhs));
	}
	template<typename T>
	Matrix operator-(T const& lhs, Matrix const& rhs)
	{
		return Matrix(rhs._cols, rhs._rows, static_cast<Complex>(lhs)) - rhs;
	}
	Matrix operator-(Matrix const& lhs, Matrix const& rhs);

	template<typename T>
	Matrix operator*(Matrix const& lhs, T const& rhs){
		Matrix result(lhs);
		for( auto& elem : result._data ){
			elem *= rhs;
		}
		return result;
	}
	template<typename T>
	Matrix operator*(T const& lhs, Matrix const& rhs){
		return rhs * lhs;
	}
	Matrix operator*(Matrix const& lhs, Matrix const& rhs);

	template<typename T>
	Matrix operator/(Matrix const& lhs, T const& rhs){
		return lhs * (1.L/rhs);
	}

    bool same_dimension(Matrix const& lhs, Matrix const& rhs);

    std::ostream& operator<<(std::ostream& out, Matrix const& matrix);
}
#endif // Feynumeric_MATRIX_HPP