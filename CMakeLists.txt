cmake_minimum_required(VERSION 3.17)
project(Feynumeric)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} C:\\fmt)
set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} "-static")
#set(FMT_LIB C:\\fmt\\cmake-build-debug)
#if (CMAKE_BUILD_TYPE EQUAL "DEBUG")
#    set(FMT_LIB C:\\fmt\\cmake-build-debug)
#elseif (CMAKE_BUILD_TYPE EQUAL "RELEASE")
#    set(FMT_LIB C:\\fmt\\cmake-build-release)
#endif()

include_directories(include)
include_directories(C:\\fmt\\include)

add_library(Feynumeric
        include/Feynumeric/units.hpp
        include/Feynumeric/table.hpp
        source/range.cpp include/Feynumeric/range.hpp
        source/particle.cpp include/Feynumeric/particle.hpp
        source/angular_momentum.cpp include/Feynumeric/angular_momentum.hpp
        source/momentum.cpp include/Feynumeric/momentum.hpp
        source/matrix.cpp include/Feynumeric/matrix.hpp
        source/kinematics.cpp include/Feynumeric/kinematics.hpp
        source/graph.cpp source/topologies.cpp include/Feynumeric/topologies.hpp source/vertex_manager.cpp include/Feynumeric/vertex_manager.hpp source/vertex.cpp include/Feynumeric/vertex.hpp source/integrate.cpp include/Feynumeric/integrate.hpp include/Feynumeric/utility.hpp source/utility.cpp source/dirac.cpp include/Feynumeric/dirac.hpp include/Feynumeric/messages.hpp source/messages.cpp source/lorentz_index.cpp include/Feynumeric/lorentz_index.hpp include/Feynumeric/function_types.hpp source/constexpr_math.cpp include/feynumeric/constexpr_math.hpp include/feynumeric/lorentz_transformation.hpp source/lorentz_transformation.cpp include/feynumeric/sum.hpp include/feynumeric/product.hpp source/feynman_graph.cpp include/feynumeric/feynman_graph.hpp source/feynman_diagram.cpp include/feynumeric/feynman_diagram.hpp include/feynumeric/constants.hpp source/topology.cpp include/feynumeric/topology.hpp include/feynumeric/format.hpp include/feynumeric/direction.hpp source/qed.cpp include/feynumeric/qed.hpp include/feynumeric/particles.hpp source/particles.cpp include/feynumeric/types.hpp source/feynman_process.cpp include/feynumeric/feynman_process.hpp source/four_vector.cpp include/feynumeric/four_vector.hpp source/three_vector.cpp include/feynumeric/three_vector.hpp include/feynumeric/core.hpp source/feynman_rule.cpp include/feynumeric/feynman_rule.hpp include/feynumeric/contract.hpp include/feynumeric/feynumeric.hpp source/particle_manager.cpp include/feynumeric/particle_manager.hpp include/feynumeric/edge_direction.hpp source/command_line_manager.cpp include/feynumeric/command_line_manager.hpp include/feynumeric/isospin.hpp source/isospin.cpp source/timer.cpp include/feynumeric/timer.hpp)

find_package(OpenMP REQUIRED)



find_library(FMT_LIB NAMES fmt HINTS "C:\\fmt\\build")
target_link_libraries(Feynumeric PRIVATE "${FMT_LIB}" OpenMP::OpenMP_CXX)
target_include_directories(Feynumeric PRIVATE include/Feynumeric)
target_compile_options(Feynumeric PRIVATE -Werror -Wall -Wextra)



#add_subdirectory(tests)
add_subdirectory(photo_production)
add_subdirectory(electron_muon)