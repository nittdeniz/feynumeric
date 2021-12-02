#ifndef FEYNCALC_FUNCTION_TYPES_HPP
#define FEYCNALC_FUNCTION_TYPES_HPP

#include <functional>
#include <memory>
#include <variant>

namespace Feyncalc
{
    class Diagram;
    using Edge_Function         = std::function<Matrix()>;
    using Vertex_Function       = std::function<Matrix(Diagram* diagram, std::vector<Edge> const& edges)>;
    using Amplitude_Function    = std::variant<Edge_Function, Vertex_Function>;
}

#endif // FEYNCALC_FUNCTION_TYPES_HPP