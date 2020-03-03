
file(GLOB_RECURSE FILENAMES *.cpp)

add_executable(run_benchmark ${FILENAMES})

target_link_libraries(run_benchmark mockturtle)
