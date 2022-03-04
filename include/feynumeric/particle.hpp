#ifndef Feynumeric_PARTICLE_HPP
#define Feynumeric_PARTICLE_HPP

#include <any>
#include <functional>
#include <memory>
#include <map>
#include <string>

#include "angular_momentum.hpp"
#include "complex.hpp"
#include "feynman_graph.hpp"
#include "matrix.hpp"
#include "momentum.hpp"
#include "types.hpp"

namespace Feynumeric
{
    class Particle
    {
    public:
        enum class Type
        {
            Majorana,
            Particle,
            AntiParticle
        };
    private:
        std::string _name;
        Type _type;
        double _mass;
        double _charge;
        Angular_Momentum _spin;
        Angular_Momentum _isospin;
        int _parity;
        std::function<double(double)> _width_function;
        double _width;
        std::map<std::string, std::any> _user_data;

    public:
        Particle(std::string&& name, Type type, double mass = 0., double width = 0., double charge = 0., double spin = 0.);
        Particle(Particle const& copy);
        Particle& operator=(Particle const& copy);

        std::string name() const;
        double mass() const;
        int charge() const;
        Angular_Momentum spin() const;
        Angular_Momentum isospin()const;

        double width() const;
        double width(double p2) const;
        void width(std::function<double(double)> f);

        void parity(int p);
        int parity() const;

        bool is_fermion() const;
        bool is_anti_fermion() const;

        std::function<Matrix(Feynman_Graph::Edge_Ptr e, Kinematics const&)> feynman_outgoing;
        std::function<Matrix(Feynman_Graph::Edge_Ptr e, Kinematics const&)> feynman_incoming;
        std::function<Matrix(Feynman_Graph::Edge_Ptr e, Kinematics const&)> feynman_virtual;

        unsigned int n_lorentz_indices() const;

        std::any user_data(std::string key) const;

        void user_data(std::string key, std::any data);

        friend std::ostream& operator<<(std::ostream&, Particle const& p);
    };

    bool is_fermion(Particle const& particle);
    bool is_anti_fermion(Particle const& particle);
    std::ostream& operator<<(std::ostream& out, Particle const& p);
}

#endif // Feynumeric_PARTICLE_HPP