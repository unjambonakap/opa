project(opa_common C)

add_library(opa_common_hdr ./dummy.c)
target_include_directories(opa_common_hdr PUBLIC ./common_base_hdr)
