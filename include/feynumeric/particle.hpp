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
#include "format.hpp"
#include "matrix.hpp"
#include "momentum.hpp"
#include "types.hpp"

namespace Feynumeric
{
    class Particle
    {
    public:
        enum class Type : uint64_t
        {
            TrueParticle = 0x01,
            AntiParticle = 0x02,
            NeutralParticle = 0x04,

            Boson           = 0x08,
	        TrueBoson       = Boson | TrueParticle,
	        AntiBoson       = Boson | AntiParticle,
	        NeutralBoson    = Boson | NeutralParticle,

	        Fermion         = 0x10,
            TrueFermion     = Fermion | TrueParticle,
            AntiFermion     = Fermion | AntiParticle,
            NeutralFermion  = Fermion | NeutralParticle,

            Lepton          = 0x20   | Fermion,
            TrueLepton      = Lepton | TrueParticle,
            AntiLepton      = Lepton | AntiParticle,

            Meson           = 0x40  | Boson,
            TrueMeson       = Meson | TrueParticle,
            AntiMeson       = Meson | AntiParticle,
            NeutralMeson    = Meson | NeutralParticle,

            Baryon          = 0x80   | Fermion,
            TrueBaryon      = Baryon | TrueParticle,
            AntiBaryon      = Baryon | AntiParticle,
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
        bool _is_group;

        double _baryon_number;
        double _lepton_number;

        Particle_Ptr _parent;

        std::map<std::string, std::any> _user_data;

    public:
    	Particle() = default;
        Particle(std::string&& name, Type type, double mass = 0., double width = 0., double charge = 0., double spin = 0.);
        Particle(Particle const& copy);
        Particle& operator=(Particle const& copy);

        std::string name() const;
        double mass() const;
        double charge() const;
        Angular_Momentum spin() const;
        Angular_Momentum isospin()const;

        void isospin(Angular_Momentum const& i);

        void copy_parameters(Particle const& other);

        void validate() const;

        double width() const;
        double width(double p2) const;
        void width(std::function<double(double)> f);

        double baryon_number() const;
        double lepton_number() const;

        void is_group(bool b);
        bool is_group() const;

        void baryon_number(double n);
        void lepton_number(double n);

        void set_type(Type const& type);

        Particle_Ptr parent() const;

        void parity(int p);
        int parity() const;

        bool is_fermion() const;
        bool is_true_fermion() const;
        bool is_anti_fermion() const;

        std::function<Matrix(Feynman_Graph::Edge_Ptr e, Kinematics const&)> feynman_outgoing;
        std::function<Matrix(Feynman_Graph::Edge_Ptr e, Kinematics const&)> feynman_incoming;
        std::function<Matrix(Feynman_Graph::Edge_Ptr e, Kinematics const&)> feynman_virtual;

        unsigned int n_lorentz_indices() const;

        std::any user_data(std::string key) const;

        template<typename T>
        T user_data(std::string key) const
        {
	        try{
		        return std::any_cast<T>(_user_data.at(key));
	        }
	        catch( std::out_of_range const& e)
	        {
		        critical_error(FORMAT("Key <{}> is not set for particle {}\n", key, _name));
	        }
        }

        void user_data(std::string key, std::any data);

        friend std::ostream& operator<<(std::ostream&, Particle const& p);
        friend class Particle_Manager;
    };

    inline uint64_t operator&(Particle::Type const& a, Particle::Type const& b)
	{
		return static_cast<uint64_t>(a) & static_cast<uint64_t>(b);
	}

	inline bool is_set(Particle::Type const& mask, Particle::Type const& flag)
	{
    	return ( mask & flag) == static_cast<uint64_t>(flag);
	}

    bool is_fermion(Particle const& particle);
    bool is_anti_fermion(Particle const& particle);
    std::ostream& operator<<(std::ostream& out, Particle const& p);
}

#endif // Feynumeric_PARTICLE_HPP