#ifndef FEYNCALC_PARTICLE_HPP
#define FEYNCALC_PARTICLE_HPP

#include <any>
#include <functional>
#include <memory>
#include <map>
#include <string>

#include "angular_momentum.hpp"
#include "momentum.hpp"
#include "complex.hpp"

namespace Feyncalc
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

        friend bool is_fermion(Particle const& particle);
        friend bool is_anti_fermion(Particle const& particle);

        std::function<Matrix()> feynman_outgoing = [](){return Matrix();};
        std::function<Matrix()> feynman_incoming = [](){return Matrix();};
        std::function<Matrix()> feynman_virtual = [](){return Matrix();};

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

#endif // FEYNCALC_PARTICLE_HPP