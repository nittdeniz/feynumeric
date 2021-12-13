#ifndef Feynumeric_PARTICLE_HPP
#define Feynumeric_PARTICLE_HPP

#include <any>
#include <functional>
#include <memory>
#include <map>
#include <string>

#include "feynumeric/angular_momentum.hpp"
#include "feynumeric/edge.hpp"
#include "feynumeric/matrix.hpp"
#include "feynumeric/momentum.hpp"
#include "feynumeric/complex.hpp"

namespace Feynumeric
{
    using std::any;
    using std::map;
    using std::string;
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
        string _name;
        Type _type;
        double _mass;
        int _charge;
        Angular_Momentum _spin;
        Angular_Momentum _isospin;
        std::function<double(double)> _width;
        map<string, any> _user_data;

    public:
        Particle(string&& name, Type type, double mass = 0, int charge = 0, Angular_Momentum spin = 0);

        string name() const;
        double mass() const;
        int charge() const;
        Angular_Momentum spin() const;
        Angular_Momentum isospin()const;

        bool is_fermion() const;
        bool is_anti_fermion() const;

        std::function<Matrix(Edge_Ptr const& e)> feynman_outgoing;
        std::function<Matrix(Edge_Ptr const& e)> feynman_incoming;
        std::function<Matrix(Edge_Ptr const& e)> feynman_virtual;

        unsigned int n_lorentz_indices() const;

        any user_data(string key) const;

        void user_data(string key, any data);

        friend std::ostream& operator<<(std::ostream&, Particle const& p);
    };

    bool is_fermion(Particle const& particle);
    bool is_anti_fermion(Particle const& particle);
    std::ostream& operator<<(std::ostream& out, Particle const& p);



    using Particle_Ptr = std::shared_ptr<Particle>;
}

#endif // Feynumeric_PARTICLE_HPP