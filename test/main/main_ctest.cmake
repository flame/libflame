
set(CTEST_MAIN_COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${PROJECT_NAME})
if(WIN32)
    set(CTEST_WORKING_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${CMAKE_BUILD_TYPE})
else()
    set(CTEST_WORKING_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
endif()

# Added test to run main test suite
add_test(NAME main_test COMMAND ${CTEST_MAIN_COMMAND} WORKING_DIRECTORY ${CTEST_WORKING_DIR})

#Example to add further tests to ctest
add_test(NAME custom_main_test_gesv_sdcz_10x10 COMMAND ${CTEST_MAIN_COMMAND} gesv sdzc 10 10 10 10 100)

#Example to add loop based tests to ctest
# Note: in forloop based test the variable should also modify the name of the test, since 2 tests cannot have same name 
foreach(FUNCTION "gesv")
    foreach(PREC "s" "d" "c" "z") 
        foreach(SIZE_N "10")
            add_test(NAME custom_main_test_${PREC}${FUNCTION}_${SIZE_N}x${SIZE_N} COMMAND ${CTEST_MAIN_COMMAND} ${FUNCTION} ${PREC} ${SIZE_N} 10 10 10 100)
        endforeach(SIZE_N)
    endforeach(PREC)
endforeach(FUNCTION)
