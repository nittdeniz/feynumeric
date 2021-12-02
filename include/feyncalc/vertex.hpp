#ifndef FEYNCALC_VERTEX_HPP
#define FEYNCALC_VERTEX_HPP

#include <vector>

namespace Feyncalc
{
    class Vertex
    {
    private:
        std::vector<std::size_t> _edge_ids;
    public:
        Vertex(std::vector<std::size_t> const& edge_ids);
    };
}
#endif // FEYNCALC_VERTEX_HPP