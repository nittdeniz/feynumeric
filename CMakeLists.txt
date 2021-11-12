cmake_minimum_required(VERSION 3.20)
project(feyncalc)

set(CMAKE_CXX_STANDARD 20)

include_directories(include)

add_library(feyncalc
        include/feyncalc/units.hpp
        include/feyncalc/table.hpp
        source/range.cpp include/feyncalc/range.hpp
        include/feyncalc/topology.hpp
        source/diagram.cpp include/feyncalc/diagram.hpp
        source/particle.cpp include/feyncalc/particle.hpp
        source/process.cpp include/feyncalc/process.hpp
        source/angular_momentum.cpp include/feyncalc/angular_momentum.hpp
        source/momentum.cpp include/feyncalc/momentum.hpp
        source/matrix.cpp include/feyncalc/matrix.hpp
)

target_include_directories(feyncalc PRIVATE include/feyncalc)
target_compile_options(feyncalc PRIVATE -Werror -Wall -Wextra)

add_subdirectory(tests)
add_subdirectory(photo_production)