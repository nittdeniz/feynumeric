/*
#ifndef FEYNUMERIC_EDGE_MANAGER_HPP
#define FEYNUMERIC_EDGE_MANAGER_HPP

#include "edge.hpp"

#include <vector>

namespace Feynumeric
{
    class Edge_Manager
    {
    private:
        std::vector<Edge> _edges;
    public:
        void add_edge();
        Edge& get(std::size_t index) const;
    };
}

#endif // FEYNUMERIC_EDGE_MANAGER_HPP*/