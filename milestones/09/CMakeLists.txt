# Defining the executable target
add_executable(09 main.cpp)

# Linking against our MD code (propagates all dependencies)
target_link_libraries(09 PUBLIC my_md_lib)
