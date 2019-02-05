set(libraries ${libraries} libplll-static)
set(libraries_shared ${libraries_shared} libplll-shared)

include_directories(${pLLL_SOURCE_DIR}/plll/include ${PROJECT_BINARY_DIR}/plll/include)
link_directories(${pLLL_BINARY_DIR}/plll)

add_executable(example-01 src/example-01.cpp)
add_executable(example-02 src/example-02.cpp)

target_link_libraries(example-01 ${libraries})
target_link_libraries(example-02 ${libraries})
