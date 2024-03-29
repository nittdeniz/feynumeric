#include <iomanip>
#include <iostream>
#include "feynumeric/format.hpp"
#include "feynumeric/matrix.hpp"
#include "feynumeric/messages.hpp"
#include "feynumeric/utility.hpp"

namespace Feynumeric
{
	Matrix::Matrix(Complex const& data)
	: _rows(1)
	, _cols(1)
	, _data({data})
	{

	}

	std::size_t Matrix::index(size_t row, std::size_t col) const
    {
        return row * _cols + col;
    }

    Matrix::Matrix(size_t rows, std::size_t cols, std::vector<Complex>&& data)
    : _rows(rows)
    , _cols(cols)
    , _data(std::move(data))
    {
        if( _data.size() != _rows * _cols )
        {
            _data.insert(_data.end(), _rows * _cols - _data.size(), Complex());
        }
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

    Matrix::Matrix(size_t rows, std::size_t cols, Complex diagonal)
    : _rows(rows)
    , _cols(cols)
    {
        if( rows != cols )
        {
        	critical_error(FORMAT("Diagonal Matrix constructors needs rows == cols (given: {} != {}})", rows, cols));
        }
        _data.resize(rows*cols);
        for( std::size_t i = 0; i < _rows; i++ )
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
        for( std::size_t i = 0; i < _data.size(); ++i )
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
        for( std::size_t i = 0; i < _data.size(); ++i )
        {
            _data[i] -= rhs._data[i];
        }
        return *this;
    }

    Matrix &Matrix::operator*=(const Matrix &rhs)
    {
        auto result = (*this) * rhs;
        _data = result._data;
        _cols = result._cols;
        _rows = result._rows;
        return *this;
    }

	Matrix& Matrix::operator*=(Complex const& rhs)
	{
		for( auto& elem : _data )
        {
            elem *= rhs;
        }
        return *this;
	}

    Matrix& Matrix::operator/=(Complex const& rhs)
    {
        for( auto& elem : _data )
        {
            elem /= rhs;
        }
        return *this;
    }

	Matrix& Matrix::operator*=(double const& rhs)
	{
		for( auto& elem : _data )
		{
			elem *= rhs;
		}
		return *this;
	}

	Complex Matrix::det() const{
		auto copy = Matrix(*this);
		double scale = 1;
//		std::cout << std::setprecision(15) << copy << "\n";
		for( std::size_t i = 0; i < copy._rows-1; ++i ){
			// chose max value as pivot point: https://en.wikipedia.org/wiki/Pivot_element
			auto norm = std::abs(copy(i, i));
			std::pair<std::size_t, double> max{i, norm};
			std::pair<std::size_t, double> min{i, norm};

			for( std::size_t j = i+1; j < copy._rows; ++j ){
				auto temp = std::abs(copy(j, i));
				if( temp > max.second ){
					max.first = j;
					max.second = temp;
				}
				if( temp < min.second ){
					min.first = j;
					min.second = temp;
				}
			}

			if( i != max.first ){
//				std::cout << "swap: " << i << " <> " << max.first << "\n";
				copy.swap_row(i, max.first);
				scale *= -1.;
			}

			if( almost_identical(copy(i, i), 0.) ){ // not full if max element is zero
				return 0.;
			}

			// rescale everything
			/*
			auto average = (min.second + max.second)/2.;
			std::cout << "average: " << average << "\n";
			copy /= average;
			scale *= int_pow(average, _rows);
			 */

//			std::cout << std::setprecision(15) << copy << "\n";

			for( std::size_t j = i+1; j < copy._rows; ++j ){
				if( almost_identical(copy(j, i), 0.) ){
//					std::cout << "almost 0: " << copy(j, i) << "\n";
					continue;
				}
				Complex factor = -copy(j, i) / copy(i, i);
				copy._data[j*copy._cols + i] = 0;
				for( std::size_t k = i+1; k < copy._cols; ++k ){
					copy._data[j*copy._cols + k] = copy._data[i*copy._cols + k] * factor + copy._data[j*copy._cols + k];
				}
			}
//			std::cout << std::setprecision(15) << copy << "\n\n\n";
		}
		Complex det = 1.;
		for( std::size_t i = 0; i < copy._rows; ++i ){
			det *= copy(i ,i);
		}
//		std::cout << "scale: " << scale << "\n";
//		std::cout << "det: " << std::fixed << std::setw(15) << std::setprecision(15) << scale * det << "\n";
		return scale * det;
	}

	void Matrix::swap_col(std::size_t i, std::vector<Complex>& col){
		if( col.size() != _rows ){
			critical_error("Dimensions for column swapping do not match.");
		}
		for( std::size_t j = 0; j < _rows; ++j ){
			std::swap(_data[j*_cols +i], col[j]);
		}
	}

	void Matrix::swap_row(std::size_t i, std::size_t j){
		if( i == j ) return;
		for( std::size_t a = 0; a < _cols; ++a ){
			std::swap(_data[i*_cols + a], _data[j*_cols + a]);
		}
	}

	Matrix& Matrix::operator*=(int const& rhs)
	{
		for( auto& elem : _data )
		{
			elem *= rhs;
		}
		return *this;
	}

	Matrix& Matrix::operator/=(double const& rhs)
	{
		for( auto& elem : _data )
		{
			elem /= rhs;
		}
		return *this;
	}

	Matrix& Matrix::operator/=(int const& rhs)
	{
		for( auto& elem : _data )
		{
			elem /= rhs;
		}
		return *this;
	}

	bool same_dimension(const Matrix &lhs, const Matrix &rhs)
    {
        return lhs._cols == rhs._cols && lhs._rows == rhs._rows;
    }

    Matrix &Matrix::apply(const std::function<Complex(Complex const&)> &f)
    {
        for( auto& elem : _data )
        {
            elem = f(elem);
        }
        return *this;
    }

    Complex &Matrix::operator()(size_t i, std::size_t j)
    {
        return _data[index(i, j)];
    }

    Complex Matrix::at(size_t i, std::size_t j) const
    {
        return _data.at(index(i, j));
    }

    Complex Matrix::try_as_complex() const
    {
        if( _rows == 1 && _cols == 1 )
        {
            return _data[0];
        }
        throw dimension_exception();
    }

    std::vector<Complex>::const_iterator Matrix::cbegin() const
    {
        return _data.cbegin();
    }

    std::vector<Complex>::const_iterator Matrix::cend() const
    {
        return _data.cend();
    }

    std::vector<Complex>::iterator Matrix::begin()
    {
        return _data.begin();
    }

    std::vector<Complex>::iterator Matrix::end()
    {
        return _data.end();
    }

    Matrix Matrix::T() const
    {
        Matrix copy(_cols, _rows);
        for( std::size_t i = 0; i < _rows; ++i )
        {
            for( std::size_t j = 0; j < _cols; ++j )
            {
                copy(j, i) = at(index(i, j));
            }
        }
        return copy;
    }

    std::size_t Matrix::elements() const
    {
		return _cols * _rows;
    }

    Complex &Matrix::operator[](size_t i)
    {
        return _data[i];
    }

    std::size_t Matrix::n_rows() const
    {
        return _rows;
    }

    std::size_t Matrix::n_cols() const
    {
        return _cols;
    }

    std::size_t Matrix::n_total() const
    {
        return _rows * _cols;
    }

    std::ostream &operator<<(std::ostream &out, const Matrix &matrix)
    {
		out << "{";
        for( std::size_t i = 0; i < matrix._rows; i++ )
        {
        	if( i > 0 ) out << ",";
            out << "{";
            for( std::size_t j = 0; j < matrix._cols; j++ )
            {

                out << ((j>0)?",":"") << matrix.at(i, j);
            }
            out << "}\n";
        }
        out << "}";
        return out;
    }

    Complex const& Matrix::at(size_t i) const
    {
        return _data[i];
    }

    Matrix operator*(const Matrix &lhs, const Matrix &rhs)
    {
        if( lhs._cols == 1 && lhs._rows == 1 )
        {
            return lhs.at(0,0) * rhs;
        }
        if( rhs._cols == 1 && rhs._rows == 1 )
        {
            return rhs.at(0,0) * lhs;
        }
        if( lhs._cols != rhs._rows )
        {
        	throw Matrix::dimension_exception();
        }
        Matrix result(lhs._rows, rhs._cols);
        result._data.resize(lhs._rows * rhs._cols);
        for( std::size_t i = 0; i < lhs._rows; i++ )
        {
            for( std::size_t j = 0; j < rhs._cols; j++ )
            {
                for( std::size_t k = 0; k < lhs._cols; k++ )
                {
                    result(i, j) += lhs.at(i, k) * rhs.at(k, j);
                }
            }
        }
        return result;
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
        	throw Matrix::dimension_exception();
        }
        Matrix result(lhs);
        std::size_t i = 0;
        for( auto& elem : result._data )
        {
            elem -= rhs._data[i++];
        }
        return result;
    }

    Matrix operator+(const Matrix &lhs, const Matrix &rhs)
    {
        Matrix result(lhs);
        std::size_t i = 0;
        for( auto& elem : result._data )
        {
            elem += rhs._data[i++];
        }
        return result;
    }

    bool operator==(Matrix const& lhs, Matrix const& rhs)
    {
        if( lhs._cols != rhs._cols || lhs._rows != rhs._rows )
        {
            return false;
        }
        bool is_same = true;
        for( std::size_t i = 0; i < lhs._data.size() && is_same; i++ )
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
