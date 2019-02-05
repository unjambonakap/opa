set(libraries ${libraries} libplll-static)
set(libraries_shared ${libraries_shared} libplll-shared)

include_directories(${pLLL_SOURCE_DIR}/plll/include ${PROJECT_BINARY_DIR}/plll/include ${pLLL_SOURCE_DIR}/plll/include-internal)
link_directories(${pLLL_BINARY_DIR}/plll)

add_executable(plll src/plll.cpp)
add_executable(pruningoptimize src/pruningoptimize.cpp)

target_link_libraries(plll ${libraries})
target_link_libraries(pruningoptimize ${libraries})

install(TARGETS plll DESTINATION bin)
