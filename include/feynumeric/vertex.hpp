#ifndef Feynumeric_VERTEX_HPP
#define Feynumeric_VERTEX_HPP

#include <vector>

namespace Feynumeric
{
    class Vertex
    {
    private:
        std::vector<std::size_t> _edge_ids;
    public:
        Vertex(std::vector<std::size_t> const& edge_ids);
    };
}
#endif // Feynumeric_VERTEX_HPP