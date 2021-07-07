add_executable(xFlame_gbench      Flame_gbench_main.cc  \
./LIN/*.cc ./EIG/*.cc ./EIG/NEIG/*.cc ./EIG/SEIG/*.cc ./EIG/ORTHO/*.cc ./EIG/COS/*.cc ./LLS/*.cc ./SVD/*cc ./AUX/*.cc  \ )

target_link_libraries(xFlame_gbench      lapack)

add_test(xFlame_gbench      ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/xFlame_gbench)
