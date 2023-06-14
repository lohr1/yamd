add_executable(04 main.cpp lj_direct_summation.cpp lj_direct_summation.h)
message(STATUS "Including header.cpp")

target_link_libraries(04 PUBLIC my_md_lib)
