add_executable(xlapacke_test_main      lapacke_test_main.cc lapacke_example_aux.c lapacke_example_aux.h \
./LIN/*.cc ./EIG/*.cc ./EIG/NEIG/*.cc ./EIG/SEIG/*.cc ./EIG/ORTHO/*.cc ./EIG/COS/*.cc ./LLS/*.cc ./SVD/*cc ./AUX/*.cc  \ )

target_link_libraries(xlapacke_test_main      lapacke)

add_test(lapacke_test_main      ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/xlapacke_test_main)
