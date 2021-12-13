#ifndef FeynumericPP_UNITS_HPP
#define FeynumericPP_UNITS_HPP

/*! \mainpage My Personal Index Page
 *
 * \section intro_sec Introduction
 *
 * This is the introduction.
 *
 * \section install_sec Installation
 *
 * \subsection step1 Step 1: Opening the box
 *
 * etc...
 */

namespace Feynumeric
{
    /**
     * @brief Energy Units for easy conversion
     */
    namespace Units
    {
        /**
         * @brief Tera Electron Volts
         * @param value in TeV
         * @return returns internal representation (in GeV)
         */
        constexpr long double operator"" _TeV(long double value)
        {
            return value * 1.e3;
        }

        constexpr long double operator"" _GeV(long double value)
        {
            return value;
        }

        constexpr long double operator"" _MeV(long double value)
        {
            return value * 1.e-3;
        }

        constexpr long double operator"" _keV(long double value)
        {
            return value * 1.e-6;
        }

        constexpr long double operator"" _eV(long double value)
        {
            return value * 1.e-9;
        }
    }
}

#endif // FeynumericPP_UNITS_HPP