#include "vertex.hpp"

namespace Feynumeric
{
    Vertex::Vertex(std::size_t id, const std::vector<Edge *> &edges)
    : _id(id)
    , _edges(edges)
    {

    }
}