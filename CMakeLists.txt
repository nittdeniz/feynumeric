cmake_minimum_required(VERSION 3.20)
project(Feynumeric)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} C:\\fmt)
#if (CMAKE_BUILD_TYPE EQUAL "DEBUG")
#    set(FMT_LIB C:\\fmt\\cmake-build-debug)
#elseif (CMAKE_BUILD_TYPE EQUAL "RELEASE")
#    set(FMT_LIB C:\\fmt\\cmake-build-release)

include_directories(include)
include_directories(C:\\fmt\\include)

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
        source/graph.cpp source/topologies.cpp include/Feynumeric/topologies.hpp source/vertex_manager.cpp include/Feynumeric/vertex_manager.hpp source/vertex.cpp include/Feynumeric/vertex.hpp source/integrate.cpp include/Feynumeric/integrate.hpp include/Feynumeric/utility.hpp source/utility.cpp source/dirac.cpp include/Feynumeric/dirac.hpp source/edge.cpp include/Feynumeric/edge.hpp include/Feynumeric/messages.hpp source/messages.cpp source/lorentz_index.cpp include/Feynumeric/lorentz_index.hpp include/Feynumeric/function_types.hpp source/constexpr_math.cpp include/feynumeric/constexpr_math.hpp include/feynumeric/lorentz_transformation.hpp source/lorentz_transformation.cpp include/feynumeric/sum.hpp include/feynumeric/product.hpp source/edge_manager.cpp include/feynumeric/edge_manager.hpp source/feynman_graph.cpp include/feynumeric/feynman_graph.hpp source/feynman_diagram.cpp include/feynumeric/feynman_diagram.hpp include/feynumeric/constants.hpp source/topology.cpp include/feynumeric/topology.hpp include/feynumeric/format.hpp include/feynumeric/direction.hpp source/qed.cpp include/feynumeric/qed.hpp include/feynumeric/particles.hpp source/particles.cpp include/feynumeric/types.hpp source/feynman_process.cpp include/feynumeric/feynman_process.hpp include/feynumeric/four_vector.cpp include/feynumeric/four_vector.hpp include/feynumeric/three_vector.cpp include/feynumeric/three_vector.hpp)

find_package(OpenMP REQUIRED)

find_library(FMT_LIB fmt)
target_link_libraries(Feynumeric PRIVATE "${FMT_LIB}" OpenMP::OpenMP_CXX)
target_include_directories(Feynumeric PRIVATE include/Feynumeric)
target_compile_options(Feynumeric PRIVATE -Werror -Wall -Wextra)



add_subdirectory(tests)
add_subdirectory(photo_production)
add_subdirectory(electron_muon)