#include "vertex_manager.hpp"
#include "messages.hpp"
#include <iostream>

namespace Feynumeric
{
    void Vertex_Manager::add_vertex( Vertex_Function const& function, std::vector<std::pair<Particle_Ptr, Vertex_Manager::Direction>> const& particles)
    {
        using std::abort;
        using std::cerr;
        using std::string;
        using std::vector;

        vector<String_Direction_Pair> pairs;
        pairs.reserve(particles.size());
        for( auto const& pair : particles )
        {
            pairs.emplace_back(pair.first->name(), pair.second);
        }

        const auto keys = generate_keys(pairs);
        _vertex_rules[keys.first][keys.second] = function;
    }

    Vertex_Manager::Vertex_Manager()
    {

    }

    Vertex_Manager::Vertex_Manager(const Vertex_Manager &other)
    : _vertex_rules(other._vertex_rules)
    {

    }

    Vertex_Manager &Vertex_Manager::operator=(const Vertex_Manager &other)
    {
        _vertex_rules = other._vertex_rules;
        return *this;
    }

    Vertex_Function Vertex_Manager::get_vertex(std::vector<String_Direction_Pair> const& pairs) const
    {
        const auto keys = generate_keys(pairs);
        try
        {
            auto particle_matches = _vertex_rules.at(keys.first);
            for( auto const& match : particle_matches )
            {
                if( match.first.size() != keys.second.size() )
                {
                    continue;
                }
                bool found = true;
                for( std::size_t i = 0; i < keys.second.size(); ++i )
                {
                    auto lhs = static_cast<unsigned int>(match.first[i]);
                    auto rhs = static_cast<unsigned int>(keys.second[i]);
                    auto result = lhs & rhs;
                    found &=  result > 0;

                }
                if( found )
                {
                    return match.second;
                }
            }
        }
        catch( std::out_of_range const& exception )
        {
            std::cerr << "Did not find vertex: " << keys.first << "\n";
            abort();
        }
        std::cerr << "Vertex does not support particle configuration.\n";
        abort();
    }

    Vertex_Function Vertex_Manager::get_vertex_function(Vertex_Id vertex_id, const vector<Edge> &edges) const
    {
        std::vector<String_Direction_Pair> pairs;
        for( auto const& edge : edges )
        {
            Direction d;
            if( edge.a() == vertex_id )
            {
                d = Direction::OUT;
            }
            else if( edge.b() == vertex_id )
            {
                d = Direction::IN;
            }
            else
            {
                critical_error("Edge " + edge.to_string() + " is not associated with vertex: " + std::to_string(static_cast<std::size_t>(vertex_id)) + ".");
            }
            pairs.emplace_back(edge.particle()->name(), d);
        }
        return get_vertex(pairs);
    }

    std::pair<std::string, std::vector<Vertex_Manager::Direction>> Vertex_Manager::generate_keys(std::vector<String_Direction_Pair> pairs) const
    {
        std::sort(pairs.begin(), pairs.end(), [](String_Direction_Pair const& lhs, String_Direction_Pair const& rhs){
            auto const& l_name = lhs.first;
            auto const& r_name = rhs.first;
            auto const& l_direction = lhs.second;
            auto const& r_direction = rhs.second;
            return l_name < r_name
                || (l_name == r_name && l_direction < r_direction);
        });

        string ID;
        vector<Direction> directions;
        directions.reserve(pairs.size());
        for( auto const& pair : pairs )
        {
            ID += pair.first;
            directions.push_back(pair.second);
        }
        return {ID, directions};
    }

}


