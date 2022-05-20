MESSAGE(STATUS "Included Tests")

set(Catch2_DIR "C:/Catch2/lib/cmake/Catch2")

#set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH} ")


message(STATUS ${CMAKE_PREFIX_PATH})

include_directories(include)
add_executable(test main.cpp)

target_link_libraries(test PRIVATE Feynumeric)
target_include_directories(test PRIVATE include ../include)

find_package(Catch2 3 REQUIRED CONFIG)
target_link_libraries(test PRIVATE Catch2::Catch2WithMain)