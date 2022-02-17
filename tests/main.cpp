#define CONFIG_CATCH_MAIN
#include <catch2/catch_test_macros.hpp>

#include <iostream>
#include <Feynumeric/angular_momentum.hpp>
#include <Feynumeric/matrix.hpp>
#include <Feynumeric/constexpr_math.hpp>
#include <Feynumeric/dirac.hpp>
#include <Feynumeric/units.hpp>
#include <Feynumeric/momentum.hpp>

TEST_CASE( "is_valid_spin", "[angular momentum]" ) {
    REQUIRE( Feynumeric::Angular_Momentum::is_valid_spin(0.) == true );
    REQUIRE( Feynumeric::Angular_Momentum::is_valid_spin(0.5) == true );
    REQUIRE( Feynumeric::Angular_Momentum::is_valid_spin(1) == true );
    REQUIRE( Feynumeric::Angular_Momentum::is_valid_spin(1.5) == true );
    REQUIRE( Feynumeric::Angular_Momentum::is_valid_spin(-1.5) == false );
    REQUIRE( Feynumeric::Angular_Momentum::is_valid_spin(2.3) == false);
    REQUIRE( Feynumeric::Angular_Momentum::is_valid_spin(-2.3) == false);
}

TEST_CASE( "is_half_odd_integer", "[angular momentum]" ) {
    Feynumeric::Angular_Momentum one_half(0.5);
    Feynumeric::Angular_Momentum one(1);
    Feynumeric::Angular_Momentum three_halves(1.5);
    REQUIRE( one_half.is_half_odd_integer() );
    REQUIRE( !one.is_half_odd_integer() );
    REQUIRE( three_halves.is_half_odd_integer() );
}

TEST_CASE( "matrix diagonal constructor", "[matrix]" ) {
    Feynumeric::Matrix a(4, 4, 3.5);
    REQUIRE( a == a);
    REQUIRE( !(a == 2*a) );
    REQUIRE( a != 2*a );
    REQUIRE( a == Feynumeric::Matrix(4, 4, {{3.5,0,0,0, 0,3.5,0,0, 0,0,3.5,0, 0,0,0,3.5}}));
}

TEST_CASE( "matrix_addition", "[matrix]" ){
    Feynumeric::Matrix a(3, 3, {1,2,3, 4,5,6, 7,8,9});
    Feynumeric::Matrix b(3, 3, {2,4,6, 8,10,12, 14,16,18});
    Feynumeric::Matrix c(3, 3, {11,2,3, 4,15,6, 7,8,19});
    REQUIRE( a + a == b );
    REQUIRE( a + 10 == c);
}

TEST_CASE( "matrix_subtraction", "[matrix]" ){
    Feynumeric::Matrix a(3, 3, {1,2,3, 4,5,6, 7,8,9});
    Feynumeric::Matrix b(3, 3, {2,4,6, 8,10,12, 14,16,18});
    Feynumeric::Matrix c(3, 3, {11,2,3, 4,15,6, 7,8,19});
    REQUIRE( a - a == Feynumeric::Matrix(3, 3, 0) );
    REQUIRE( b - a == a);
    REQUIRE( -(10-c) == a);
}

TEST_CASE( "matrix_operator()", "[matrix]"){
    Feynumeric::Matrix a(3, 3, {1,2,3, 4,5,6, 7,8,9});
    Feynumeric::Matrix b(3, 3, {1,2,13, 4,5,6, 7,8,9});
    REQUIRE( a != b );
    a(0, 2) += 10;
    REQUIRE( a == b );
}

TEST_CASE( "matrix_multiplication", "[matrix]"){
    Feynumeric::Matrix a(1, 4, {1,2,3,4});
    Feynumeric::Matrix b(4, 1, {5,6,7,8});

    REQUIRE( a*b == Feynumeric::Matrix(1, 1, std::vector<Feynumeric::Complex>{70}));
    REQUIRE( b*a == Feynumeric::Matrix(4, 4, {5,10,15,20, 6,12,18,24, 7,14,21,28, 8,16,24,32}));
}

TEST_CASE( "matrix_division", "[matrix]" ){
    Feynumeric::Matrix a(3, 3, {1,2,3, 4,5,6, 7,8,9});
    Feynumeric::Matrix b(3, 3, {2,4,6, 8,10,12, 14,16,18});

    REQUIRE(b/2. == a);
}

TEST_CASE( "matrix transposition", "[matrix]"){
    Feynumeric::Matrix a(3, 3, {1,2,3, 4,5,6, 7,8,9});
    Feynumeric::Matrix b(3, 3, {1,4,7, 2,5,8, 3,6,9});

//    REQUIRE(a.T() == b);
}

TEST_CASE("clebsch_gordan", "[math]"){
    using namespace Feynumeric;
    REQUIRE(is_almost_equal(clebsch_gordan(0,0,0,0,0,0), 1.));
    REQUIRE(is_almost_equal(clebsch_gordan(2,0,3,0,5,0), constexpr_sqrt(10./21.)));
    REQUIRE(is_almost_equal(clebsch_gordan(2,1,3,0,5,1), 2*constexpr_sqrt(2./21.)));
    REQUIRE(is_almost_equal(clebsch_gordan(2,0,3,1,5,1), constexpr_sqrt(3./7.)));
    REQUIRE(is_almost_equal(clebsch_gordan(2,1,3,1,5,2), constexpr_sqrt(1./2.)));
    REQUIRE(is_almost_equal(clebsch_gordan(0.5, 0.5, 0.5, -0.5, 1, 0), constexpr_sqrt(1./2.)));
    REQUIRE(is_almost_equal(clebsch_gordan(0.5, 0.5, 0.5, -0.5, 0, 0), constexpr_sqrt(1./2.)));
    REQUIRE(is_almost_equal(clebsch_gordan(2.5,0.5,1.5,0.5,2,1), -5./(2*constexpr_sqrt(21.))));
    // exceptions
//    REQUIRE([](){
//        try{ clebsch_gordan(-1,0,0,0,0,0);return false; }
//        catch( std::domain_error& e ){ return true; }}()
//    );
//    REQUIRE([](){
//        try{ clebsch_gordan(1,0,1,0,3,0);return false; }
//        catch( std::domain_error& e ){ return true; }}()
//    );
//    REQUIRE([](){
//        try{ clebsch_gordan(1,0,2,0,3,1);return false; }
//        catch( std::domain_error& e ){ return true; }}()
//    );
}

TEST_CASE("gamma^2 == 4", "[dirac]"){
	using namespace Feynumeric;
	REQUIRE( GA[0] * GAC[0] + GA[1] * GAC[1] + GA[2] * GAC[2] + GA[3] * GAC[3] == Matrix(4,4,4));
}

TEST_CASE("Dirac Spinors Completeness", "[dirac]"){
	using namespace Feynumeric;
	using namespace Feynumeric::Units;
	Particle_Ptr Muon_Minus   = std::make_shared<Particle>("Muon_-", Particle::Type::Particle, 105.6583745_MeV, -1, 0.5);
	Four_Momentum p(0.3_GeV, Muon_Minus->mass(), std::cos(0.3), std::cos(0.2));
	Angular_Momentum_Ptr s1 = std::make_shared<Angular_Momentum>(0.5, 0.5);
	Angular_Momentum_Ptr s2 = std::make_shared<Angular_Momentum>(0.5, -0.5);
	auto lhs = u(Muon_Minus, p, s1, {}) * ubar(Muon_Minus, p, s1, {}) + u(Muon_Minus, p, s2, {}) * ubar(Muon_Minus, p, s2, {});
	auto rhs = GS(p) + Muon_Minus->mass();
	for( std::size_t i = 0; i < 16; ++i )
	{
		REQUIRE( almost_identical(lhs.at(i), rhs.at(i), 0.1) );
	}
}

TEST_CASE("Spin 1 Polarisation Vectors Completeness", "[dirac]"){
	using namespace Feynumeric;

	Four_Momentum p(3, 4, 0.1234, 0.5918273);

	Lorentz_Index_Ptr mu = std::make_shared<Lorentz_Index>();
	Lorentz_Index_Ptr nu = std::make_shared<Lorentz_Index>();

	Angular_Momentum_Ptr s1p = std::make_shared<Angular_Momentum>(1, 1);
	Angular_Momentum_Ptr s0 = std::make_shared<Angular_Momentum>(1, 0);
	Angular_Momentum_Ptr s1m = std::make_shared<Angular_Momentum>(1, -1);

	Matrix result(4, 4);
	Matrix compare(4, 4);
	for( int i = 0; i < 4; ++i )
	{
		for( int j = 0; j < 4; ++j )
		{
			auto temp1 = epsilon(p, s1p, {mu}) * epsilon_star(p, s1p, {nu});
			auto temp2 = epsilon(p, s0, {mu}) * epsilon_star(p, s0, {nu});
			auto temp3 = epsilon(p, s1m, {mu}) * epsilon_star(p, s1m, {nu});
			result(i, j) = (temp1 + temp2 + temp3).try_as_complex();
			auto temp4 = -MT[*mu][*nu];
			auto temp5 = FV(p, mu) * FV(p, nu);
			auto temp6 = temp5 / p.dot();
			compare(i, j) =  temp4 + temp6;
			++(*nu);
		}
		++(*mu);
	}
	for( std::size_t i = 0; i < 16; ++i )
	{
		REQUIRE( almost_identical(result.at(i), compare.at(i)) );
	}

}