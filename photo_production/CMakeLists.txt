add_executable(coupling_constants
        coupling_constants.cpp
        source/effective_lagrangian_model.cpp include/particles.hpp

        ../electron_muon/include/vertices.hpp ../electron_muon/source/vertices.cpp include/effective_lagrangian_model.hpp include/form_factors.hpp source/form_factors.cpp)

add_executable(playground playground.cpp)

add_executable(photo_production
        cross_sections.cpp
        source/effective_lagrangian_model.cpp include/particles.hpp

        ../electron_muon/include/vertices.hpp ../electron_muon/source/vertices.cpp include/effective_lagrangian_model.hpp
        include/form_factors.hpp source/form_factors.cpp)

target_include_directories(photo_production PRIVATE include)
target_include_directories(coupling_constants PRIVATE include)
target_include_directories(playground PRIVATE include)

target_link_libraries(photo_production PRIVATE Feynumeric)
target_link_libraries(coupling_constants PRIVATE Feynumeric)
target_link_libraries(playground PRIVATE Feynumeric)