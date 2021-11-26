#include "vertex_manager.hpp"
#include <iostream>

namespace Feyncalc
{
    void Vertex_Manager::add_vertex( Vertex_Manager::Vertex_Function const& function, std::vector<std::pair<Particle_Ptr, Vertex_Manager::Direction>> const& particles)
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

    Vertex_Manager::Vertex_Function Vertex_Manager::get_vertex(std::vector<String_Direction_Pair> const& pairs) const
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
                    found &= static_cast<unsigned int>(match.first[i]) &
                             static_cast<unsigned int>(keys.second[i]);
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

    std::pair<std::string, std::vector<Vertex_Manager::Direction>> Vertex_Manager::generate_keys(std::vector<String_Direction_Pair> pairs) const
    {
        std::sort(pairs.begin(), pairs.end(), [](String_Direction_Pair const& lhs, String_Direction_Pair const& rhs){
            return lhs.first < rhs.first;
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


