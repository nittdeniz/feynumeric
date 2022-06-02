add_executable(coupling_constants
        coupling_constants.cpp
        source/effective_lagrangian_model.cpp include/particles.hpp

        ../electron_muon/include/vertices.hpp ../electron_muon/source/vertices.cpp include/effective_lagrangian_model.hpp include/form_factors.hpp source/form_factors.cpp include/write_map.hpp source/write_map.cpp include/coupling_constants.hpp)

add_executable(playground playground.cpp source/effective_lagrangian_model.cpp source/form_factors.cpp)

add_executable(dyson_factors dyson_factors.cpp source/effective_lagrangian_model.cpp source/form_factors.cpp source/write_map.cpp)

add_executable(cross_sections
        cross_sections.cpp
        source/effective_lagrangian_model.cpp include/particles.hpp

        ../electron_muon/include/vertices.hpp ../electron_muon/source/vertices.cpp include/effective_lagrangian_model.hpp
        include/form_factors.hpp source/form_factors.cpp)

target_include_directories(cross_sections PRIVATE include)
target_include_directories(coupling_constants PRIVATE include)
target_include_directories(playground PRIVATE include)
target_include_directories(dyson_factors PRIVATE include)

target_link_libraries(cross_sections PRIVATE Feynumeric)
target_link_libraries(coupling_constants PRIVATE Feynumeric)
target_link_libraries(playground PRIVATE Feynumeric)
target_link_libraries(dyson_factors PRIVATE Feynumeric OpenMP::OpenMP_CXX)