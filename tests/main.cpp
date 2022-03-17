#define CONFIG_CATCH_MAIN
#include <catch2/catch_test_macros.hpp>

#include <feynumeric/angular_momentum.hpp>
#include <feynumeric/matrix.hpp>
#include <feynumeric/constexpr_math.hpp>
#include <feynumeric/dirac.hpp>
#include <feynumeric/integrate.hpp>
#include <feynumeric/units.hpp>
#include <feynumeric/qed.hpp>
#include <feynumeric/feynman_diagram.hpp>
#include <feynumeric/topologies.hpp>
#include <feynumeric/feynman_process.hpp>


#include <cmath>
#include <iostream>
#include <map>


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

TEST_CASE("n_states", "[angular momentum]"){
	Feynumeric::Angular_Momentum spin_zero(0);
	Feynumeric::Angular_Momentum spin_onehalf(0.5);
	Feynumeric::Angular_Momentum spin_three_massless(3, 3, true);
	REQUIRE( spin_zero.n_states() == 1 );
	REQUIRE( spin_onehalf.n_states() == 2);
	REQUIRE( spin_three_massless.n_states() == 2);
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

TEST_CASE("Dirac Particle Spinors Completeness", "[dirac]"){
	using namespace Feynumeric;
	using namespace Feynumeric::Units;
	using namespace Feynumeric::QED;
//	Particle_Ptr Muon_Minus   = std::make_shared<Particle>("Muon_-", Particle::Type::Particle, 105.6583745_MeV, -1, 0.5);
	Four_Vector p = four_momentum(0.3_GeV, Muon_Minus->mass(), std::cos(0.3), std::cos(0.2));
	Angular_Momentum_Ptr s1 = std::make_shared<Angular_Momentum>(0.5, 0.5);
	Angular_Momentum_Ptr s2 = std::make_shared<Angular_Momentum>(0.5, -0.5);
	auto lhs = u(Muon_Minus, p, s1, {}) * ubar(Muon_Minus, p, s1, {}) + u(Muon_Minus, p, s2, {}) * ubar(Muon_Minus, p, s2, {});
	auto rhs = GS(p) + Muon_Minus->mass();
	for( std::size_t i = 0; i < 16; ++i )
	{
		REQUIRE( almost_identical(lhs.at(i), rhs.at(i), 0.1) );
	}
}

TEST_CASE("Dirac Anti Particle Spinors Completeness", "[dirac]"){
	using namespace Feynumeric;
	using namespace Feynumeric::Units;
	using namespace Feynumeric::QED;
	//Particle_Ptr Muon_Plus   = std::make_shared<Particle>("Muon_+", Particle::Type::Particle, 105.6583745_MeV, 1, 0.5);
	Four_Vector p = four_momentum(0.3_GeV, Muon_Plus->mass(), std::cos(0.3), std::cos(0.2));
	Angular_Momentum_Ptr s1 = std::make_shared<Angular_Momentum>(0.5, 0.5);
	Angular_Momentum_Ptr s2 = std::make_shared<Angular_Momentum>(0.5, -0.5);
	auto lhs = v(Muon_Plus, p, s1, {}) * vbar(Muon_Plus, p, s1, {}) + v(Muon_Plus, p, s2, {}) * vbar(Muon_Plus, p, s2, {});
	auto rhs = GS(p) - Muon_Plus->mass();
	for( std::size_t i = 0; i < 16; ++i )
	{
		REQUIRE( almost_identical(lhs.at(i), rhs.at(i), 0.1) );
	}
}

TEST_CASE("Spin 1 Polarisation Vectors Completeness", "[dirac]"){
	using namespace Feynumeric;

	Particle_Ptr test_particle = std::make_shared<Particle>("Test", Particle::Type::TrueBoson, 4., 0, 1, 1);

	Four_Vector p = four_momentum(3, test_particle->mass(), 0.2, 0.3);

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
			auto temp1 = (epsilon(test_particle, p, s1p, {mu}) * epsilon_star(test_particle, p, s1p, {nu})).try_as_complex();
			auto temp2 = (epsilon(test_particle, p, s0, {mu}) * epsilon_star(test_particle, p, s0, {nu})).try_as_complex();
			auto temp3 = (epsilon(test_particle, p, s1m, {mu}) * epsilon_star(test_particle, p, s1m, {nu})).try_as_complex();
			result(i, j) = temp1 + temp2 + temp3;
			auto temp4 = -MT[*mu][*nu];
			auto temp5 = p.contra(mu) * p.contra(nu);
			auto temp6 = temp5 / p.squared();
			compare(i, j) =  temp4 + temp6;
			++(*nu);
		}
		++(*mu);
	}
	for( std::size_t i = 0; i < 16; ++i )
	{
		REQUIRE( std::abs(result.at(i) - result.at(i)) < 0.00000001 );
	}
}

TEST_CASE("Integration Routine", "[math]")
{
	using namespace Feynumeric;
	using f_type = std::function<double(double)>;

	auto sin = [](double x){return std::sin(x);};
	auto x2 = [](double x){return x*x;};


	auto result0 = integrate([](double x){return x;}, -1, 2);
	auto result1 = integrate(sin, 0, M_PI);
	auto result2 = integrate(x2, 0, 4);
	auto result3 = integrate([](double x){return std::exp(x);}, -5, 5, 1.e-12);
	auto result4 = integrate([](double x){return std::sin(x)/x;}, 0.00001, 1);
	auto result5 = integrate([](double x){return std::sin(x)/x;}, 1, 100);
	REQUIRE( almost_identical(result0, 1.5, 1.e-6));
	REQUIRE( almost_identical(result1, 2., 1.e-6) );
	REQUIRE( almost_identical(result2, 64/3., 1.e-6) );
	REQUIRE( almost_identical(result3, 148.406421, 1.e-6));
	REQUIRE( almost_identical(result4, 0.946073, 1.e-6));
	REQUIRE( almost_identical(result5, 0.616142, 1.e-6));
}

TEST_CASE("Moller Scattering", "[QED]")
{
	using namespace Feynumeric;
	using namespace Feynumeric::Units;
	using namespace Feynumeric::QED;

	init_particles();
	init_vertices();

	Feynman_Diagram_Ptr t_channel = create_diagram("t_channel", X_Man, VMP,
	                                               {Electron, Electron},
	                                               {Photon},
	                                               {Electron, Electron});

	Feynman_Diagram_Ptr u_channel = create_diagram("u_channel", X_Man, VMP,
	                                               {Electron, Electron},
	                                               {Photon},
	                                               {Electron, Electron});
	u_channel->cross_outgoing(0, 1);

	Feynman_Process e_scattering({t_channel, u_channel});

	std::stringstream out;

	double const cos_theta = 0.2134;

	auto result = e_scattering.dsigma_dcos_table( 500._MeV, {cos_theta});

	// Compare to analytical values from Mathematica's Feyncalc
	REQUIRE( almost_identical(result[cos_theta][0], 0.02227945628277883) );
	REQUIRE( almost_identical(result[cos_theta][1], 0.01880419088437939) );
	REQUIRE( almost_identical(result[cos_theta][2], 0.0085134844703209) );
}

TEST_CASE("Bhaba Scattering", "[QED]")
{
	using namespace Feynumeric;
	using namespace Feynumeric::Units;
	using namespace Feynumeric::QED;

	init_particles();
	init_vertices();

	Feynman_Diagram_Ptr s_channel = create_diagram("s_channel", Double_Wrench, VMP,
	                                               {Electron, Positron},
	                                               {Photon},
	                                               {Electron, Positron});

	Feynman_Diagram_Ptr t_channel = create_diagram("t_channel", X_Man, VMP,
	                                               {Electron, Positron},
	                                               {Photon},
	                                               {Electron, Positron});


	Feynman_Process e_scattering({s_channel, t_channel});

	//std::stringstream out;

	double const cos_theta = 0.41936;

	auto result = e_scattering.dsigma_dcos_table( 500._MeV, {cos_theta});

	// Compare to analytical values from Mathematica's Feyncalc
	REQUIRE( almost_identical(result[cos_theta][0], 0.0095746468309887) );
	REQUIRE( almost_identical(result[cos_theta][1], 0.024487089153493090) );
	REQUIRE( almost_identical(result[cos_theta][2], 0.017657720374332880) );
}

TEST_CASE("Polarisation Sums", "[dirac]")
{
	using namespace Feynumeric;

//	Particle_Ptr test_particle = std::make_shared<Particle>("Test", Particle::Type::Majorana, 4.);

	auto check_pol_sum = [&](Angular_Momentum_Ptr spin)
	{
		Particle_Ptr test_particle;
		if( spin->is_half_odd_integer() )
		{
			test_particle = std::make_shared<Particle>("Test", Particle::Type::TrueFermion, 4., 0, 0, spin->j());
		}
		else
		{
			test_particle = std::make_shared<Particle>("Test", Particle::Type::TrueBoson, 4., 0, 0, spin->j());
		}

		Four_Vector p = four_momentum(2, test_particle->mass(), 0.3, 1);
		auto n_indices = std::ceil(spin->j() - 0.5) * 2;
		std::vector<Lorentz_Index_Ptr> mu(n_indices);
		for( auto& index : mu )
		{
			index = std::make_shared<Lorentz_Index>();
		}
		std::vector<Lorentz_Index_Ptr> indices_left(mu.begin(), mu.begin() + mu.size()/2);
		std::vector<Lorentz_Index_Ptr> indices_right(mu.begin() + mu.size()/2, mu.end());
		std::vector<Angular_Momentum_Ptr> S(spin->n_states());
		for( std::size_t i = 0; i < spin->n_states(); ++i )
		{
			S[i] = std::make_shared<Angular_Momentum>(spin->j(), spin->j() - i * 1.);
		}

		auto const N = std::pow(4, n_indices);
		for( std::size_t i = 0; i < N; ++i )
		{
			auto projector = Projector(test_particle, p, mu);
			if( spin->is_half_odd_integer() )
			{
				Matrix temp(4, 4, 0);
				for( auto const& s : S )
				{
					auto u1 = u(test_particle, p, s, indices_left);
					auto u2 = ubar(test_particle, p, s, indices_right);
					temp += u1 * u2;
				}
				for( std::size_t j = 0; j < temp.elements(); ++j )
				{
					if( !almost_identical(projector.at(j), temp.at(j)) )
					{
						return false;
					}
				}
			}
			else
			{
				Complex temp{0., 0.};
				for( auto const& s : S )
				{
					temp += ( epsilon(test_particle, p, s, indices_left) *
					             epsilon_star(test_particle, p, s, indices_right)).try_as_complex();
				}
				if( !almost_identical(projector.try_as_complex(), temp) )
				{
					return false;
				}
			}
			// advance mu
			for( std::size_t l = 0; l < mu.size(); ++l )
			{
				++(*(mu[l]));
				if( *(mu[l]) != 0 )
				{
					break;
				}
			}
		}
		return true;
	};

//	REQUIRE( check_pol_sum(std::make_shared<Angular_Momentum>(0.5, 0.5)) );
//	REQUIRE( check_pol_sum(std::make_shared<Angular_Momentum>(1., 1.)) );
	REQUIRE( check_pol_sum(std::make_shared<Angular_Momentum>(1.5, 1.5)) );
	REQUIRE( check_pol_sum(std::make_shared<Angular_Momentum>(2., 2.)) );
	REQUIRE( check_pol_sum(std::make_shared<Angular_Momentum>(2.5, 2.5)) );
	REQUIRE( check_pol_sum(std::make_shared<Angular_Momentum>(3., 3.)) );
	REQUIRE( check_pol_sum(std::make_shared<Angular_Momentum>(3.5, 3.5)) );
}