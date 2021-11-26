#include "dirac.hpp"

namespace Feyncalc
{
    [[maybe_unused]] Matrix GS(const Matrix &a)
    {
        return Matrix(a);
    }

    [[maybe_unused]] Matrix dirac_sigma(const Matrix &a, const Matrix &b)
    {
        return Complex(0, 1)/2. * (a*b - b*a);
    }

    [[maybe_unused]] Matrix epsilon1(const Matrix &p, const Angular_Momentum &lambda)
    {
        switch( static_cast<int>(lambda.m()) )
        {
            case 1:
            case 0:
            case -1:
            default:
                break;
        }
        return Matrix(p);
    }

    [[maybe_unused]] Matrix epsilon2(const Matrix &p, const Angular_Momentum &lambda)
    {
        switch( static_cast<int>(lambda.m()) )
        {
            case 2:
            case 1:
            case 0:
            case -1:
            case -2:
            default:
                break;
        }
        return Matrix(p);
    }

    [[maybe_unused]] Matrix u12(const Matrix &p, const Angular_Momentum &s)
    {
        switch( static_cast<int>(2*s.m()) )
        {
            case 1:
            case -1:
            default:
                break;
        }
        return Matrix(p);
    }

    [[maybe_unused]] Matrix ubar12(const Matrix &p, const Angular_Momentum &s)
    {
        switch( static_cast<int>(2*s.m()) )
        {
            case 1:
            case -1:
            default:
                break;
        }
        return Matrix(p);
    }

    [[maybe_unused]] Matrix u32(const Matrix &p, const Angular_Momentum &s)
    {
        switch( static_cast<int>(2*s.m()) )
        {
            case 3:
            case 1:
            case -1:
            case -3:
            default:
                break;
        }
        return Matrix(p);
    }

    [[maybe_unused]] Matrix ubar32(const Matrix &p, const Angular_Momentum &s)
    {
        switch( static_cast<int>(2*s.m()) )
        {
            case 3:
            case 1:
            case -1:
            case -3:
            default:
                break;
        }
        return Matrix(p);
    }

    [[maybe_unused]] Matrix u52(const Matrix &p, const Angular_Momentum &s)
    {
        switch( static_cast<int>(2*s.m()) )
        {
            case 5:
            case 3:
            case 1:
            case -1:
            case -3:
            case -5:
            default:
                break;
        }
        return Matrix(p);
    }

    [[maybe_unused]] Matrix ubar52(const Matrix &p, const Angular_Momentum &s)
    {
        switch( static_cast<int>(2*s.m()) )
        {
            case 5:
            case 3:
            case 1:
            case -1:
            case -3:
            case -5:
            default:
                break;
        }
        return Matrix(p);
    }

    [[maybe_unused]] Matrix v12(const Matrix &p, const Angular_Momentum &s)
    {
        switch( static_cast<int>(2*s.m()) )
        {
            case 1:
            case -1:
            default:
                break;
        }
        return Matrix(p);
    }

    [[maybe_unused]] Matrix vbar12(const Matrix &p, const Angular_Momentum &s)
    {
        switch( static_cast<int>(2*s.m()) )
        {
            case 1:
            case -1:
            default:
                break;
        }
        return Matrix(p);
    }

    [[maybe_unused]] Matrix v32(const Matrix &p, const Angular_Momentum &s)
    {
        switch( static_cast<int>(2*s.m()) )
        {
            case 3:
            case 1:
            case -1:
            case -3:
            default:
                break;
        }
        return Matrix(p);
    }

    [[maybe_unused]] Matrix vbar32(const Matrix &p, const Angular_Momentum &s)
    {
        switch( static_cast<int>(2*s.m()) )
        {
            case 3:
            case 1:
            case -1:
            case -3:
            default:
                break;
        }
        return Matrix(p);
    }

    [[maybe_unused]] Matrix v52(const Matrix &p, const Angular_Momentum &s)
    {
        switch( static_cast<int>(2*s.m()) )
        {
            case 5:
            case 3:
            case 1:
            case -1:
            case -3:
            case -5:
            default:
                break;
        }
        return Matrix(p);
    }

    [[maybe_unused]] Matrix vbar52(const Matrix &p, const Angular_Momentum &s)
    {
        switch( static_cast<int>(2*s.m()) )
        {
            case 5:
            case 3:
            case 1:
            case -1:
            case -3:
            case -5:
            default:
                break;
        }
        return Matrix(p);
    }
}