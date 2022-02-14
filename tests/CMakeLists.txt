MESSAGE(STATUS "Included Tests")

include_directories(include)
add_executable(test
        main.cpp
        source/particles.cpp
        source/test_angular_momentum.cpp)

target_link_libraries(test PRIVATE Feynumeric)
target_include_directories(test PRIVATE include ../include)


find_package(Catch2 3 REQUIRED)
target_link_libraries(test PRIVATE Catch2::Catch2WithMain)