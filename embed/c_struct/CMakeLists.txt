project(c_struct C)


set(SRC ./ring_buffer.c)
add_library(c_struct STATIC ${SRC})
target_include_directories(c_struct PUBLIC ./inc)
target_link_libraries(c_struct PUBLIC ${C_STRUCT_DEPS})


