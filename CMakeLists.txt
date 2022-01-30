cmake_minimum_required(VERSION 3.20)
project(Feynumeric)

set(CMAKE_CXX_STANDARD 20)

include_directories(include)

add_library(Feynumeric
        include/Feynumeric/units.hpp
        include/Feynumeric/table.hpp
        source/range.cpp include/Feynumeric/range.hpp
        include/Feynumeric/graph.hpp
        source/diagram.cpp include/Feynumeric/diagram.hpp
        source/particle.cpp include/Feynumeric/particle.hpp
        source/process.cpp include/Feynumeric/process.hpp
        source/angular_momentum.cpp include/Feynumeric/angular_momentum.hpp
        source/momentum.cpp include/Feynumeric/momentum.hpp
        source/matrix.cpp include/Feynumeric/matrix.hpp
        source/kinematics.cpp include/Feynumeric/kinematics.hpp
        source/graph.cpp source/topologies.cpp include/Feynumeric/topologies.hpp source/vertex_manager.cpp include/Feynumeric/vertex_manager.hpp source/vertex.cpp include/Feynumeric/vertex.hpp source/integrate.cpp include/Feynumeric/integrate.hpp include/Feynumeric/utility.hpp source/utility.cpp source/dirac.cpp include/Feynumeric/dirac.hpp source/edge.cpp include/Feynumeric/edge.hpp include/Feynumeric/messages.hpp source/messages.cpp source/lorentz_index.cpp include/Feynumeric/lorentz_index.hpp include/Feynumeric/function_types.hpp include/feynumeric/constexpr_math.cpp include/feynumeric/constexpr_math.hpp include/feynumeric/boost.hpp source/boost.cpp include/feynumeric/sum.hpp include/feynumeric/product.hpp)

target_include_directories(Feynumeric PRIVATE include/Feynumeric)
target_compile_options(Feynumeric PRIVATE -Werror -Wall -Wextra)

add_subdirectory(tests)
add_subdirectory(photo_production)
add_subdirectory(electron_muon)