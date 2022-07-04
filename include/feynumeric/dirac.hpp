#ifndef Feynumeric_DIRAC_HPP
#define Feynumeric_DIRAC_HPP


#include "lorentz_transformation.hpp"
#include "kinematics.hpp"
#include "matrix.hpp"
#include "momentum.hpp"

#include <array>
#include <functional>
#include <vector>

namespace Feynumeric
{
	class Angular_Momentum;
	class Graph_Edge;
	class Lorentz_Index;
	class Particle;
	
    extern std::array<Matrix, 4> GA;
    extern std::array<Matrix, 4> GAC;
    extern Matrix GA5;
    extern std::array<std::array<double, 4>, 4> MT;


    Matrix GS(Four_Vector const& p);
	Matrix GS(Matrix const& matrix);

	Matrix dirac_sigma(std::shared_ptr<Lorentz_Index> mu, std::shared_ptr<Lorentz_Index> nu);
	Matrix dirac_sigma(std::shared_ptr<Lorentz_Index> mu, Four_Vector const& p);
	Matrix dirac_sigma(Four_Vector const& p, std::shared_ptr<Lorentz_Index> nu);

	Matrix dirac_sigmac(std::shared_ptr<Lorentz_Index> mu, std::shared_ptr<Lorentz_Index> nu);
	Matrix dirac_sigmac(std::shared_ptr<Lorentz_Index> mu, Four_Vector const& p);
	Matrix dirac_sigmac(Four_Vector const& p, std::shared_ptr<Lorentz_Index> nu);




	Matrix O32(Four_Vector const& p, std::shared_ptr<Lorentz_Index> const& mu, std::shared_ptr<Lorentz_Index> const& nu);
	Matrix O32(Four_Vector const& p, Four_Vector const& q, std::shared_ptr<Lorentz_Index> const& nu);
	Matrix O32(Four_Vector const& p, std::shared_ptr<Lorentz_Index> const& mu, Four_Vector const& q);
	Matrix O32c(Four_Vector const& p, std::shared_ptr<Lorentz_Index> const& mu, std::shared_ptr<Lorentz_Index> const& nu);
	Matrix O32c(Four_Vector const& p, Four_Vector const& q, std::shared_ptr<Lorentz_Index> const& nu);
	Matrix O32c(Four_Vector const& p, std::shared_ptr<Lorentz_Index> const& mu, Four_Vector const& q);

    Matrix O(Four_Vector const& p, std::vector<Lorentz_Index_Ptr> const& mu, std::vector<Lorentz_Index_Ptr> lambda);
    Matrix O(Four_Vector const& p, std::vector<Four_Vector> const& Ps, std::vector<Lorentz_Index_Ptr> lambda);

    Matrix Oc(Four_Vector const& p, std::vector<Four_Vector> const& Ps, std::vector<Lorentz_Index_Ptr> lambda);
    Matrix Oc(Four_Vector const& p, std::vector<Lorentz_Index_Ptr> const& mu, std::vector<Lorentz_Index_Ptr> lambda);

    Matrix u(std::shared_ptr<Graph_Edge> edge_ptr, Kinematics const& kin);
    Matrix u(std::shared_ptr<Particle> const& P, Four_Vector const& p, std::shared_ptr<Angular_Momentum> s, std::vector<std::shared_ptr<Lorentz_Index>> const& lorentz_indices);
    Matrix ubar(std::shared_ptr<Graph_Edge> edge_ptr, Kinematics const& kin);
    Matrix ubar(std::shared_ptr<Particle> const& P, Four_Vector const& p, std::shared_ptr<Angular_Momentum> s, std::vector<std::shared_ptr<Lorentz_Index>> const& lorentz_indices);

	Matrix v(std::shared_ptr<Graph_Edge> edge_ptr, Kinematics const& kin);
	Matrix v(std::shared_ptr<Particle> const& P, Four_Vector const& p, std::shared_ptr<Angular_Momentum> s, std::vector<std::shared_ptr<Lorentz_Index>> const& lorentz_indices);
	Matrix vbar(std::shared_ptr<Graph_Edge> edge_ptr, Kinematics const& kin);
	Matrix vbar(std::shared_ptr<Particle> const& P, Four_Vector const& p, std::shared_ptr<Angular_Momentum> s, std::vector<std::shared_ptr<Lorentz_Index>> const& lorentz_indices);

	Matrix spinor(std::shared_ptr<Particle> const& P, Four_Vector const& p, std::shared_ptr<Angular_Momentum> s, std::vector<std::shared_ptr<Lorentz_Index>> const& lorentz_indices);

	Complex epsilon(std::shared_ptr<Angular_Momentum> const& spin, double q, double mass, double cos_theta, double cos_phi, std::vector<std::shared_ptr<Lorentz_Index>> const& mus);
    Matrix epsilon(std::shared_ptr<Graph_Edge> edge_ptr, Kinematics const& kin);
    Matrix epsilon_star(std::shared_ptr<Graph_Edge> edge_ptr, Kinematics const& kin);

    Matrix epsilon(std::shared_ptr<Particle> const& P, Four_Vector const& p, std::shared_ptr<Angular_Momentum> s, std::vector<std::shared_ptr<Lorentz_Index>> const& lorentz_indices);
    Matrix epsilon_star(std::shared_ptr<Particle> const& P, Four_Vector const& p, std::shared_ptr<Angular_Momentum> s, std::vector<std::shared_ptr<Lorentz_Index>> const& lorentz_indices);

    Matrix Projector(std::shared_ptr<Graph_Edge> edge_ptr, Kinematics const& kin);
    Matrix Propagator(std::shared_ptr<Graph_Edge> edge_ptr, Kinematics const& kin);

    Matrix Projector(std::shared_ptr<Particle> const& P, const Four_Vector &p, const std::vector<std::shared_ptr<Lorentz_Index>> &lorentz_indices);
	Matrix Projector(std::shared_ptr<Particle> const& P, Angular_Momentum const& spin, const Four_Vector &p, const std::vector<std::shared_ptr<Lorentz_Index>> &lorentz_indices);
    Matrix Propagator(std::shared_ptr<Particle> const& P, const Four_Vector &p, const std::vector<std::shared_ptr<Lorentz_Index>> &lorentz_indices);
}

#endif // Feynumeric_DIRAC_HPP


