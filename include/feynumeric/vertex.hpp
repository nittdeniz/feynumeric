#ifndef Feynumeric_VERTEX_HPP
#define Feynumeric_VERTEX_HPP

#include <vector>

namespace Feynumeric
{
    class Edge;
    class Vertex
    {
    private:
        std::size_t _id;
        std::vector<Edge*> _edges;
    public:
        Vertex(std::size_t id, std::vector<Edge*> const& edges);
    };
}
#endif // Feynumeric_VERTEX_HPP