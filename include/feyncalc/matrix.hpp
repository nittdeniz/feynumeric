#ifndef FEYNCALC_MATRIX_HPP
#define FEYNCALC_MATRIX_HPP

#include "complex.hpp"
#include <vector>

namespace Feyncalc
{
    using std::size_t;
    using std::vector;
    class Matrix
    {
    private:
        size_t _rows, _cols;
        vector<Complex> _data;
        size_t index(size_t row, size_t col) const;
    public:
        class dimension_exception : public std::exception
        {
        public:
            dimension_exception(){}
        };

        Matrix(size_t rows = 0, size_t cols = 0, vector<Complex>&& data = {});
        Matrix(size_t rows, size_t cols, Complex diagonal);
        Matrix(Matrix const& other);
        Matrix(Matrix&& other);
        Matrix& operator=(Matrix const& other);
        Matrix& operator=(Matrix&& other);

        Complex& operator()(size_t i, size_t j);
        Complex at(size_t i, size_t j) const;

        Matrix& operator+=(Matrix const& rhs);
        Matrix& operator-=(Matrix const& rhs);
        Matrix& operator*=(Matrix const& rhs);

        Complex try_as_complex() const;

        friend Matrix operator+(Matrix const& lhs, Matrix const& rhs);
        friend Matrix operator+(Matrix const& lhs, int rhs);
        friend Matrix operator+(int lhs, Matrix const& rhs);
        friend Matrix operator+(Matrix const& lhs, double rhs);
        friend Matrix operator+(double lhs, Matrix const& rhs);
        friend Matrix operator+(Matrix const& lhs, Complex rhs);
        friend Matrix operator+(Complex lhs, Matrix const& rhs);

        friend Matrix operator-(Matrix const& lhs);

        friend Matrix operator-(Matrix const& lhs, Matrix const& rhs);
        friend Matrix operator-(Matrix const& lhs, int rhs);
        friend Matrix operator-(int lhs, Matrix const& rhs);
        friend Matrix operator-(Matrix const& lhs, double rhs);
        friend Matrix operator-(double lhs, Matrix const& rhs);
        friend Matrix operator-(Matrix const& lhs, Complex rhs);
        friend Matrix operator-(Complex lhs, Matrix const& rhs);

        friend Matrix operator*(Matrix const& lhs, Matrix const& rhs);

        friend Matrix operator*(Matrix const& lhs, int rhs);
        friend Matrix operator*(int lhs, Matrix const& rhs);
        friend Matrix operator*(Matrix const& lhs, double rhs);
        friend Matrix operator*(double lhs, Matrix const& rhs);
        friend Matrix operator*(Matrix const& lhs, Complex rhs);
        friend Matrix operator*(Complex lhs, Matrix const& rhs);

        friend Matrix operator/(Matrix const& lhs, int rhs);
        friend Matrix operator/(Matrix const& lhs, double rhs);
        friend Matrix operator/(Matrix const& lhs, Complex rhs);

        friend bool operator==(Matrix const& lhs, Matrix const& rhs);
        friend bool operator!=(Matrix const& lhs, Matrix const& rhs);

        friend bool same_dimension(Matrix const& lhs, Matrix const& rhs);
        friend std::ostream& operator<<(std::ostream& out, Matrix const& matrix);
    };

    Matrix operator+(Matrix const& lhs, Matrix const& rhs);
    Matrix operator+(Matrix const& lhs, int rhs);
    Matrix operator+(int lhs, Matrix const& rhs);
    Matrix operator+(Matrix const& lhs, double rhs);
    Matrix operator+(double lhs, Matrix const& rhs);
    Matrix operator+(Matrix const& lhs, Complex rhs);
    Matrix operator+(Complex lhs, Matrix const& rhs);

    Matrix operator-(Matrix const& lhs);

    Matrix operator-(Matrix const& lhs, Matrix const& rhs);
    Matrix operator-(Matrix const& lhs, int rhs);
    Matrix operator-(int lhs, Matrix const& rhs);
    Matrix operator-(Matrix const& lhs, double rhs);
    Matrix operator-(double lhs, Matrix const& rhs);
    Matrix operator-(Matrix const& lhs, Complex rhs);
    Matrix operator-(Complex lhs, Matrix const& rhs);

    Matrix operator*(Matrix const& lhs, Matrix const& rhs);

    Matrix operator*(Matrix const& lhs, int rhs);
    Matrix operator*(int lhs, Matrix const& rhs);
    Matrix operator*(Matrix const& lhs, double rhs);
    Matrix operator*(double lhs, Matrix const& rhs);
    Matrix operator*(Matrix const& lhs, Complex rhs);
    Matrix operator*(Complex lhs, Matrix const& rhs);

    Matrix operator/(Matrix const& lhs, int rhs);
    Matrix operator/(Matrix const& lhs, double rhs);
    Matrix operator/(Matrix const& lhs, Complex rhs);

    bool same_dimension(Matrix const& lhs, Matrix const& rhs);

    std::ostream& operator<<(std::ostream& out, Matrix const& matrix);
}
#endif // FEYNCALC_MATRIX_HPP