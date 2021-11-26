#include <iostream>
#include "feyncalc/matrix.hpp"

namespace Feyncalc
{
    using std::cerr;
    size_t Matrix::index(size_t row, size_t col) const
    {
        return row * _cols + col;
    }

    Matrix::Matrix(size_t rows, size_t cols, vector<Complex>&& data)
    : _rows(rows)
    , _cols(cols)
    , _data(std::move(data))
    {

    }

    Matrix::Matrix(const Matrix &other)
    : _rows(other._rows)
    , _cols(other._cols)
    , _data(other._data)
    {
    }

    Matrix::Matrix(Matrix &&other)
    : _rows(std::move(other._rows))
    , _cols(std::move(other._cols))
    , _data(std::move(other._data))
    {

    }

    Matrix::Matrix(size_t rows, size_t cols, Complex diagonal)
    : _rows(rows)
    , _cols(cols)
    {
        if( rows != cols )
        {
            cerr << "Diagonal Matrix constructors needs rows == cols (given: " << rows << " != " << cols << ")\n";
            abort();
        }
        _data.resize(rows*cols);
        for( size_t i = 0; i < _rows; i++ )
        {
            operator()(i, i) = diagonal;
        }
    }

    Matrix &Matrix::operator=(const Matrix &other)
    {
        _rows = other._rows;
        _cols = other._cols;
        _data = other._data;
        return *this;
    }

    Matrix &Matrix::operator=(Matrix &&other)
    {
        _rows = std::move(other._rows);
        _cols = std::move(other._cols);
        _data = std::move(other._data);
        return *this;
    }

    Matrix &Matrix::operator+=(const Matrix &rhs)
    {
        if( !same_dimension(*this, rhs) )
        {
            throw dimension_exception();
        }
        for( size_t i = 0; i < _data.size(); ++i )
        {
            _data[i] += rhs._data[i];
        }
        return *this;
    }

    Matrix &Matrix::operator-=(const Matrix &rhs)
    {
        if( !same_dimension(*this, rhs) )
        {
            throw dimension_exception();
        }
        for( size_t i = 0; i < _data.size(); ++i )
        {
            _data[i] -= rhs._data[i];
        }
        return *this;
    }

    Matrix &Matrix::operator*=(const Matrix &rhs)
    {
        *this = (*this) * rhs;
        return *this;
    }


    bool same_dimension(const Matrix &lhs, const Matrix &rhs)
    {
        return lhs._cols == rhs._cols && lhs._rows == rhs._rows;
    }


    Complex &Matrix::operator()(size_t i, size_t j)
    {
        return _data[index(i, j)];
    }

    Complex Matrix::at(size_t i, size_t j) const
    {
        return _data.at(index(i, j));
    }

    Complex Matrix::try_as_complex() const
    {
        if( _rows == 1 && _cols == 1 )
        {
            return _data[0];
        }
        cerr << "rows: " << _rows << "\tcols: " << _cols << "\n";
        throw dimension_exception();
    }

    std::ostream &operator<<(std::ostream &out, const Matrix &matrix)
    {
        for( size_t i = 0; i < matrix._rows; i++ )
        {
            out << "(";
            for( size_t j = 0; j < matrix._cols; j++ )
            {

                out << ((j>0)?"\t":"") << matrix.at(i, j);
            }
            out << ")\n";
        }
        return out;
    }

    Matrix operator*(const Matrix &lhs, const Matrix &rhs)
    {
        if( lhs._cols != rhs._rows )
        {
            cerr << "lhs._cols != rhs._rows @operator*(Matrix, matrix)\n";
            cerr << "LHS:\n";
            cerr << lhs << "\n";
            cerr << "RHS:\n";
            cerr << rhs << "\n";
            abort();
        }
        Matrix result(lhs._rows, rhs._cols);
        result._data.resize(lhs._rows * rhs._cols);
        for( size_t i = 0; i < lhs._rows; i++ )
        {
            for( size_t j = 0; j < rhs._cols; j++ )
            {
                for( size_t k = 0; k < lhs._cols; k++ )
                {
                    result(i, j) += lhs.at(i, k) * rhs.at(k, j);
                }
            }
        }
        return result;
    }

    Matrix operator*(const Matrix &lhs, int rhs)
    {
        Matrix result(lhs);
        for( auto& elem : result._data )
        {
            elem *= static_cast<double>(rhs);
        }
        return result;
    }

    Matrix operator*(int lhs, const Matrix &rhs)
    {
        return operator*(rhs, lhs);
    }

    Matrix operator*(const Matrix &lhs, double rhs)
    {
        Matrix result(lhs);
        for( auto& elem : result._data )
        {
            elem *= rhs;
        }
        return result;
    }

    Matrix operator*(double lhs, const Matrix &rhs)
    {
        return operator*(rhs, lhs);
    }

    Matrix operator*(const Matrix &lhs, Complex rhs)
    {
        Matrix result(lhs);
        for( auto& elem : result._data )
        {
            elem *= rhs;
        }
        return result;
    }

    Matrix operator*(Complex lhs, const Matrix &rhs)
    {
        return operator*(rhs, lhs);
    }

    Matrix operator-(const Matrix &lhs)
    {
        Matrix result(lhs);
        for( auto& elem : result._data )
        {
            elem = -elem;
        }
        return result;
    }

    Matrix operator-(const Matrix &lhs, const Matrix &rhs)
    {
        if( !same_dimension(lhs, rhs) )
        {
            cerr << "Dimensions of Matrix operator- do not match.\n";
            cerr << "lhs:\n";
            cerr << lhs << "\n";
            cerr << "rhs:\n";
            cerr << rhs << "\n";
            abort();
        }
        Matrix result(lhs);
        size_t i = 0;
        for( auto& elem : result._data )
        {
            elem -= rhs._data[i++];
        }
        return result;
    }

    Matrix operator-(const Matrix &lhs, int rhs)
    {
        return operator+(lhs, -rhs);
    }

    Matrix operator-(int lhs, const Matrix &rhs)
    {
        return operator+(-rhs, lhs);
    }

    Matrix operator-(const Matrix &lhs, double rhs)
    {
        return operator+(lhs, -rhs);
    }

    Matrix operator-(double lhs, const Matrix &rhs)
    {
        return operator+(-rhs, lhs);
    }

    Matrix operator-(const Matrix &lhs, Complex rhs)
    {
        return operator+(lhs, -rhs);
    }

    Matrix operator-(Complex lhs, const Matrix &rhs)
    {
        return operator+(-rhs, lhs);
    }

    Matrix operator/(const Matrix &lhs, int rhs)
    {
        return operator*(lhs, 1./rhs);
    }

    Matrix operator/(const Matrix &lhs, double rhs)
    {
        return operator*(lhs, 1./rhs);
    }

    Matrix operator/(const Matrix &lhs, Complex rhs)
    {
        return operator*(lhs, 1./rhs);
    }

    Matrix operator+(const Matrix &lhs, const Matrix &rhs)
    {
        Matrix result(lhs);
        size_t i = 0;
        for( auto& elem : result._data )
        {
            elem += rhs._data[i++];
        }
        return result;
    }

    Matrix operator+(const Matrix &lhs, int rhs)
    {
        return lhs + Matrix(lhs._cols, lhs._rows, static_cast<Complex>(rhs));
    }

    Matrix operator+(int lhs, const Matrix &rhs)
    {
        return rhs + lhs;
    }

    Matrix operator+(const Matrix &lhs, double rhs)
    {
        return lhs + Matrix(lhs._cols, lhs._rows, static_cast<Complex>(rhs));
    }

    Matrix operator+(double lhs, const Matrix &rhs)
    {
        return rhs + lhs;
    }

    Matrix operator+(const Matrix &lhs, Complex rhs)
    {
        return lhs + Matrix(lhs._rows, lhs._cols, rhs);
    }

    Matrix operator+(Complex lhs, const Matrix &rhs)
    {
        return rhs + lhs;
    }

    bool operator==(Matrix const& lhs, Matrix const& rhs)
    {
        if( lhs._cols != rhs._cols || lhs._rows != rhs._rows )
        {
            return false;
        }
        bool is_same = true;
        for( size_t i = 0; i < lhs._data.size() && is_same; i++ )
        {
            is_same &= lhs._data[i] == rhs._data[i];
        }
        return is_same;
    }

    bool operator!=(Matrix const& lhs, Matrix const& rhs)
    {
        return !(lhs==rhs);
    }
}
