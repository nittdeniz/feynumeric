add_executable(photo_production
        main.cpp
        source/particles.cpp include/particles.hpp

)

target_include_directories(photo_production PRIVATE include)

target_link_libraries(photo_production PRIVATE feyncalc)