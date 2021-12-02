add_executable(photo_production
        main.cpp
        source/particles.cpp include/particles.hpp

        ../electron_muon/include/vertices.hpp ../electron_muon/source/vertices.cpp)

target_include_directories(photo_production PRIVATE include)

target_link_libraries(photo_production PRIVATE feyncalc)