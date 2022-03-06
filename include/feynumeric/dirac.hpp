#ifndef Feynumeric_DIRAC_HPP
#define Feynumeric_DIRAC_HPP

#include <array>
#include <functional>
#include <vector>

#include "angular_momentum.hpp"
#include "lorentz_transformation.hpp"
#include "feynman_graph.hpp"
#include "kinematics.hpp"
#include "matrix.hpp"
#include "momentum.hpp"
#include "particle.hpp"


namespace Feynumeric
{
    extern std::array<Matrix, 4> GA;
    extern std::array<Matrix, 4> GAC;
    extern Matrix GA5;
    extern std::array<std::array<double, 4>, 4> MT;


    Matrix GS(Four_Vector const& p);
	Matrix GS(Matrix const& matrix);

	Matrix dirac_sigma(Lorentz_Index_Ptr mu, Lorentz_Index_Ptr nu);
	Matrix dirac_sigma(Lorentz_Index_Ptr mu, Four_Vector const& p);
	Matrix dirac_sigma(Four_Vector const& p, Lorentz_Index_Ptr nu);

	Matrix dirac_sigmac(Lorentz_Index_Ptr mu, Lorentz_Index_Ptr nu);
	Matrix dirac_sigmac(Lorentz_Index_Ptr mu, Four_Vector const& p);
	Matrix dirac_sigmac(Four_Vector const& p, Lorentz_Index_Ptr nu);


//	Complex FV(Matrix const& a, Lorentz_Index_Ptr mu);
//	Complex FV(Four_Vector const& p, Lorentz_Index_Ptr mu);
//	Complex FVC(Matrix const& a, Lorentz_Index_Ptr mu);
//	Complex FVC(Four_Vector const& p, Lorentz_Index_Ptr nu);


    Matrix u(Feynman_Graph::Edge_Ptr edge_ptr, Kinematics const& kin);
    Matrix u(Particle_Ptr const& P, Four_Vector const& p, Angular_Momentum_Ptr s, std::vector<Lorentz_Index_Ptr> const& lorentz_indices);
    Matrix ubar(Feynman_Graph::Edge_Ptr edge_ptr, Kinematics const& kin);
    Matrix ubar(Particle_Ptr const& P, Four_Vector const& p, Angular_Momentum_Ptr s, std::vector<Lorentz_Index_Ptr> const& lorentz_indices);

	Matrix v(Feynman_Graph::Edge_Ptr edge_ptr, Kinematics const& kin);
	Matrix v(Particle_Ptr const& P, Four_Vector const& p, Angular_Momentum_Ptr s, std::vector<Lorentz_Index_Ptr> const& lorentz_indices);
	Matrix vbar(Feynman_Graph::Edge_Ptr edge_ptr, Kinematics const& kin);
	Matrix vbar(Particle_Ptr const& P, Four_Vector const& p, Angular_Momentum_Ptr s, std::vector<Lorentz_Index_Ptr> const& lorentz_indices);

	Matrix spinor(Particle_Ptr const& P, Four_Vector const& p, Angular_Momentum_Ptr s, std::vector<Lorentz_Index_Ptr> const& lorentz_indices);

    Matrix epsilon(Feynman_Graph::Edge_Ptr edge_ptr, Kinematics const& kin);
    Matrix epsilon_star(Feynman_Graph::Edge_Ptr edge_ptr, Kinematics const& kin);

    Matrix epsilon(Particle_Ptr const& P, Four_Vector const& p, Angular_Momentum_Ptr s, std::vector<Lorentz_Index_Ptr> const& lorentz_indices);
    Matrix epsilon_star(Particle_Ptr const& P, Four_Vector const& p, Angular_Momentum_Ptr s, std::vector<Lorentz_Index_Ptr> const& lorentz_indices);

    Matrix Projector(Feynman_Graph::Edge_Ptr edge_ptr, Kinematics const& kin, bool ignore_momentum = false);
    Matrix Propagator(Feynman_Graph::Edge_Ptr edge_ptr, Kinematics const& kin, bool ignore_momentum = false);

    Matrix Projector(Particle_Ptr const& P, const Four_Vector &p, const std::vector<Lorentz_Index_Ptr> &lorentz_indices, bool ignore_momentum = false);
	Matrix Projector(Particle_Ptr const& P, Angular_Momentum const& spin, const Four_Vector &p, const std::vector<Lorentz_Index_Ptr> &lorentz_indices, bool ignore_momentum = false);
    Matrix Propagator(Particle_Ptr const& P, const Four_Vector &p, const std::vector<Lorentz_Index_Ptr> &lorentz_indices, bool ignore_momentum = false);
}

#endif // Feynumeric_DIRAC_HPP


