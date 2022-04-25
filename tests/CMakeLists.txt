MESSAGE(STATUS "Included Tests")

set(Catch2_DIR "C:/Catch2/build")

include_directories(include)
add_executable(test main.cpp)

target_link_libraries(test PRIVATE Feynumeric)
target_include_directories(test PRIVATE include ../include)


find_package(Catch2 3 REQUIRED)
target_link_libraries(test PRIVATE Catch2::Catch2WithMain)