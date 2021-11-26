cmake_minimum_required(VERSION 3.20)
project(feyncalc)

set(CMAKE_CXX_STANDARD 20)

include_directories(include)

add_library(feyncalc
        include/feyncalc/units.hpp
        include/feyncalc/table.hpp
        source/range.cpp include/feyncalc/range.hpp
        include/feyncalc/graph.hpp
        source/diagram.cpp include/feyncalc/diagram.hpp
        source/particle.cpp include/feyncalc/particle.hpp
        source/process.cpp include/feyncalc/process.hpp
        source/angular_momentum.cpp include/feyncalc/angular_momentum.hpp
        source/momentum.cpp include/feyncalc/momentum.hpp
        source/matrix.cpp include/feyncalc/matrix.hpp
        source/kinematics.cpp include/feyncalc/kinematics.hpp
        source/graph.cpp source/topologies.cpp include/feyncalc/topologies.hpp source/vertex_manager.cpp include/feyncalc/vertex_manager.hpp source/vertex.cpp include/feyncalc/vertex.hpp source/integrate.cpp include/feyncalc/integrate.hpp include/feyncalc/utility.hpp source/utility.cpp source/dirac.cpp include/feyncalc/dirac.hpp source/edge.cpp include/feyncalc/edge.hpp)

target_include_directories(feyncalc PRIVATE include/feyncalc)
target_compile_options(feyncalc PRIVATE -Werror -Wall -Wextra)

add_subdirectory(tests)
add_subdirectory(photo_production)
add_subdirectory(electron_muon)