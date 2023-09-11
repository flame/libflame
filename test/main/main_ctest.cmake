
set(CTEST_MAIN_COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${PROJECT_NAME})
if(WIN32)
    set(CTEST_WORKING_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${CMAKE_BUILD_TYPE})
else()
    set(CTEST_WORKING_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
endif()

# Added test to run main test suite
foreach(CONFIG_TYPE "long" "medium" "short" "micro")
    add_test(NAME main_test_${CONFIG_TYPE} COMMAND ${CTEST_MAIN_COMMAND}  --config-dir=${CONFIG_TYPE} WORKING_DIRECTORY ${CTEST_WORKING_DIR})
endforeach()

#Example to add further tests to ctest
add_test(NAME custom_main_test_gesv_sdcz_10x10 COMMAND ${CTEST_MAIN_COMMAND} gesv sdzc 10 10 10 10 100)
#Performance tests for DGESVD
add_test(NAME DGESVD_SML_OPT00 COMMAND ${CTEST_MAIN_COMMAND} gesvd d N N 3 1 3 3 1 -1 1000)
set_property(TEST DGESVD_SML_OPT00 PROPERTY ENVIRONMENT "OMP_NUM_THREADS=64")
add_test(NAME DGESVD_SML_OPT01 COMMAND ${CTEST_MAIN_COMMAND} gesvd d N N 3 2 3 3 2 -1 1000)
set_property(TEST DGESVD_SML_OPT01 PROPERTY ENVIRONMENT "OMP_NUM_THREADS=64")
add_test(NAME DGESVD_SML_OPT02 COMMAND ${CTEST_MAIN_COMMAND} gesvd d N N 3 3 3 3 3 -1 1000)
set_property(TEST DGESVD_SML_OPT02 PROPERTY ENVIRONMENT "OMP_NUM_THREADS=64")
add_test(NAME DGESVD_SML_OPT03 COMMAND ${CTEST_MAIN_COMMAND} gesvd d N N 3 4 3 3 4 -1 1000)
set_property(TEST DGESVD_SML_OPT03 PROPERTY ENVIRONMENT "OMP_NUM_THREADS=64")
add_test(NAME DGESVD_SML_OPT04 COMMAND ${CTEST_MAIN_COMMAND} gesvd d N N 3 8 3 3 8 -1 1000)
set_property(TEST DGESVD_SML_OPT04 PROPERTY ENVIRONMENT "OMP_NUM_THREADS=64")
add_test(NAME DGESVD_SML_OPT05 COMMAND ${CTEST_MAIN_COMMAND} gesvd d N N 3 9 3 3 9 -1 1000)
set_property(TEST DGESVD_SML_OPT05 PROPERTY ENVIRONMENT "OMP_NUM_THREADS=64")
add_test(NAME DGESVD_SML_OPT06 COMMAND ${CTEST_MAIN_COMMAND} gesvd d N N 3 16 3 3 16 -1 100)
set_property(TEST DGESVD_SML_OPT06 PROPERTY ENVIRONMENT "OMP_NUM_THREADS=64")
add_test(NAME DGESVD_SML_OPT07 COMMAND ${CTEST_MAIN_COMMAND} gesvd d N N 3 17 3 3 17 -1 100)
set_property(TEST DGESVD_SML_OPT07 PROPERTY ENVIRONMENT "OMP_NUM_THREADS=64")
add_test(NAME DGESVD_SML_OPT08 COMMAND ${CTEST_MAIN_COMMAND} gesvd d N N 3 21 3 3 21 -1 100)
set_property(TEST DGESVD_SML_OPT08 PROPERTY ENVIRONMENT "OMP_NUM_THREADS=64")
add_test(NAME DGESVD_SML_OPT09 COMMAND ${CTEST_MAIN_COMMAND} gesvd d N N 3 23 3 3 23 -1 100)
set_property(TEST DGESVD_SML_OPT09 PROPERTY ENVIRONMENT "OMP_NUM_THREADS=64")
add_test(NAME DGESVD_SML_OPT10 COMMAND ${CTEST_MAIN_COMMAND} gesvd d N N 3 34 3 3 34 -1 100)
set_property(TEST DGESVD_SML_OPT10 PROPERTY ENVIRONMENT "OMP_NUM_THREADS=64")
add_test(NAME DGESVD_SML_OPT11 COMMAND ${CTEST_MAIN_COMMAND} gesvd d N N 3 40 3 3 40 -1 100)
set_property(TEST DGESVD_SML_OPT11 PROPERTY ENVIRONMENT "OMP_NUM_THREADS=64")
add_test(NAME DGESVD_SML_OPT12 COMMAND ${CTEST_MAIN_COMMAND} gesvd d N N 3 51 3 3 51 -1 100)
set_property(TEST DGESVD_SML_OPT12 PROPERTY ENVIRONMENT "OMP_NUM_THREADS=64")
add_test(NAME DGESVD_SML_OPT13 COMMAND ${CTEST_MAIN_COMMAND} gesvd d N N 3 53 3 3 53 -1 100)
set_property(TEST DGESVD_SML_OPT13 PROPERTY ENVIRONMENT "OMP_NUM_THREADS=64")
add_test(NAME DGESVD_SML_OPT14 COMMAND ${CTEST_MAIN_COMMAND} gesvd d N N 3 67 3 3 67 -1 100)
set_property(TEST DGESVD_SML_OPT14 PROPERTY ENVIRONMENT "OMP_NUM_THREADS=64")
add_test(NAME DGESVD_SML_OPT15 COMMAND ${CTEST_MAIN_COMMAND} gesvd d N N 3 84 3 3 84 -1 100)
set_property(TEST DGESVD_SML_OPT15 PROPERTY ENVIRONMENT "OMP_NUM_THREADS=64")
add_test(NAME DGESVD_SML_OPT16 COMMAND ${CTEST_MAIN_COMMAND} gesvd d N N 3 268 3 3 268 -1 100)
set_property(TEST DGESVD_SML_OPT16 PROPERTY ENVIRONMENT "OMP_NUM_THREADS=64")
add_test(NAME DGESVD_SML_OPT17 COMMAND ${CTEST_MAIN_COMMAND} gesvd d N N 3 364 3 3 364 -1 100)
set_property(TEST DGESVD_SML_OPT17 PROPERTY ENVIRONMENT "OMP_NUM_THREADS=64")
add_test(NAME DGESVD_SML_OPT18 COMMAND ${CTEST_MAIN_COMMAND} gesvd d S S 1 2 1 1 2 -1 1000)
set_property(TEST DGESVD_SML_OPT18 PROPERTY ENVIRONMENT "OMP_NUM_THREADS=64")
add_test(NAME DGESVD_SML_OPT19 COMMAND ${CTEST_MAIN_COMMAND} gesvd d S S 2 3 2 2 3 -1 1000)
set_property(TEST DGESVD_SML_OPT19 PROPERTY ENVIRONMENT "OMP_NUM_THREADS=64")
add_test(NAME DGESVD_SML_OPT20 COMMAND ${CTEST_MAIN_COMMAND} gesvd d S S 4 3 4 4 3 -1 1000)
set_property(TEST DGESVD_SML_OPT20 PROPERTY ENVIRONMENT "OMP_NUM_THREADS=64")
add_test(NAME DGESVD_SML_OPT21 COMMAND ${CTEST_MAIN_COMMAND} gesvd d S S 4 5 4 4 5 -1 1000)
set_property(TEST DGESVD_SML_OPT21 PROPERTY ENVIRONMENT "OMP_NUM_THREADS=64")
add_test(NAME DGESVD_SML_OPT22 COMMAND ${CTEST_MAIN_COMMAND} gesvd d S S 6 3 6 6 3 -1 1000)
set_property(TEST DGESVD_SML_OPT22 PROPERTY ENVIRONMENT "OMP_NUM_THREADS=64")
add_test(NAME DGESVD_SML_OPT23 COMMAND ${CTEST_MAIN_COMMAND} gesvd d S S 6 4 6 6 4 -1 1000)
set_property(TEST DGESVD_SML_OPT23 PROPERTY ENVIRONMENT "OMP_NUM_THREADS=64")

#Example to add loop based tests to ctest
# Note: in forloop based test the variable should also modify the name of the test, since 2 tests cannot have same name 
foreach(FUNCTION "gesv")
    foreach(PREC "s" "d" "c" "z") 
        foreach(SIZE_N "10")
            add_test(NAME custom_main_test_${PREC}${FUNCTION}_${SIZE_N}x${SIZE_N} COMMAND ${CTEST_MAIN_COMMAND} ${FUNCTION} ${PREC} ${SIZE_N} 10 10 10 100)
        endforeach(SIZE_N)
    endforeach(PREC)
endforeach(FUNCTION)
