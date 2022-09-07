add_executable(coupling_constants coupling_constants.cpp source/effective_lagrangian_model.cpp include/particles.hpp  ../electron_muon/include/vertices.hpp ../electron_muon/source/vertices.cpp include/effective_lagrangian_model.hpp include/form_factors.hpp source/form_factors.cpp include/write_map.hpp source/write_map.cpp source/couplings.cpp include/couplings.hpp)

add_executable(playground playground.cpp source/effective_lagrangian_model.cpp source/form_factors.cpp source/couplings.cpp include/couplings.hpp)

add_executable(decay_widths decay_widths.cpp source/effective_lagrangian_model.cpp source/form_factors.cpp source/write_map.cpp source/couplings.cpp include/couplings.hpp)

add_executable(dyson_factors dyson_factors.cpp source/effective_lagrangian_model.cpp source/form_factors.cpp source/write_map.cpp source/couplings.cpp include/couplings.hpp)

add_executable(cross_sections
        cross_sections.cpp
        source/effective_lagrangian_model.cpp include/particles.hpp

        ../electron_muon/include/vertices.hpp ../electron_muon/source/vertices.cpp include/effective_lagrangian_model.hpp
        include/form_factors.hpp source/form_factors.cpp source/couplings.cpp include/couplings.hpp)

add_executable(amplitude amplitude.cpp source/couplings.cpp include/couplings.hpp include/effective_lagrangian_model.hpp source/effective_lagrangian_model.cpp source/form_factors.cpp include/form_factors.hpp)

add_executable(gradient_descent gradient_descent.cpp source/couplings.cpp include/couplings.hpp include/effective_lagrangian_model.hpp source/effective_lagrangian_model.cpp source/form_factors.cpp include/form_factors.hpp)
add_executable(gradient_descent_differential gradient_descent_differential.cpp source/couplings.cpp include/couplings.hpp include/effective_lagrangian_model.hpp source/effective_lagrangian_model.cpp source/form_factors.cpp include/form_factors.hpp)

add_executable(amplitude2D amplitude2D.cpp source/couplings.cpp include/couplings.hpp include/effective_lagrangian_model.hpp source/effective_lagrangian_model.cpp include/form_factors.hpp source/form_factors.cpp)

target_include_directories(amplitude PRIVATE include)
target_include_directories(cross_sections PRIVATE include)
target_include_directories(coupling_constants PRIVATE include)
target_include_directories(playground PRIVATE include)
target_include_directories(dyson_factors PRIVATE include)
target_include_directories(decay_widths PRIVATE include)
target_include_directories(gradient_descent PRIVATE include)
target_include_directories(gradient_descent_differential PRIVATE include)
target_include_directories(amplitude2D PRIVATE include)

target_link_libraries(amplitude PRIVATE Feynumeric OpenMP::OpenMP_CXX)
target_link_libraries(cross_sections PRIVATE Feynumeric)
target_link_libraries(coupling_constants PRIVATE Feynumeric)
target_link_libraries(playground PRIVATE Feynumeric)
target_link_libraries(dyson_factors PRIVATE Feynumeric OpenMP::OpenMP_CXX)
target_link_libraries(decay_widths PRIVATE Feynumeric)
target_link_libraries(gradient_descent PRIVATE Feynumeric)
target_link_libraries(gradient_descent_differential PRIVATE Feynumeric)
target_link_libraries(amplitude2D PRIVATE Feynumeric OpenMP::OpenMP_CXX)