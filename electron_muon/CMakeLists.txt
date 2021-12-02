SET(EXECUTABLE_NAME electron_muon)

add_executable(${EXECUTABLE_NAME}
        main.cpp
        source/particles.cpp include/particles.hpp
        source/vertices.cpp include/vertices.hpp

)

target_include_directories(${EXECUTABLE_NAME} PRIVATE include)

target_link_libraries(${EXECUTABLE_NAME} PRIVATE feyncalc)