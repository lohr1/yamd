add_executable(04 main.cpp)
message(STATUS "Including header.cpp")

target_link_libraries(04 PUBLIC my_md_lib)
