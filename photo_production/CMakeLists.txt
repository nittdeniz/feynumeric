add_executable(photo_production
        main.cpp
        source/effective_lagrangian_model.cpp include/particles.hpp

        ../electron_muon/include/vertices.hpp ../electron_muon/source/vertices.cpp include/effective_lagrangian_model.hpp include/form_factors.hpp)

target_include_directories(photo_production PRIVATE include)

target_link_libraries(photo_production PRIVATE Feynumeric)