cmake_minimum_required(VERSION 3.12)
project(libflame_testsuite C)

# Force C99 standard compliance
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# ---------------------------------------------------------------------------
# Find BLAS
# ---------------------------------------------------------------------------
# Find an underlying BLAS engine (e.g., OpenBLAS, MKL, or reference BLAS)
find_package(BLAS REQUIRED)

# ---------------------------------------------------------------------------
# Find libflame
# ---------------------------------------------------------------------------
# List of environment variables to check, in order of preference
set(ENV_VAR_CANDIDATES LIBFLAME_ROOT LIBFLAME_DIR LIBFLAME_PREFIX LIBFLAME_PATH)
set(DETERMINED_ROOT "")

# Loop through candidates to find the first one that is defined
foreach(VAR_NAME ${ENV_VAR_CANDIDATES})
    if(DEFINED ENV{${VAR_NAME}} AND NOT "$ENV{${VAR_NAME}}" STREQUAL "")
        set(DETERMINED_ROOT "$ENV{${VAR_NAME}}")
        message(STATUS "Found libflame root path via environment variable: ${VAR_NAME}=${DETERMINED_ROOT}")
        break() # Stop searching once a valid path is found
    endif()
endforeach()

# Safety Check: If after the loop DETERMINED_ROOT is still empty, halt configuration
if("${DETERMINED_ROOT}" STREQUAL "")
    message(FATAL_ERROR
        "Could not locate libflame! You must export one of the following environment variables "
        "pointing to your libflame installation directory:\n"
        "  -> LIBFLAME_ROOT\n"
        "  -> LIBFLAME_DIR\n"
        "  -> LIBFLAME_PREFIX\n"
        "  -> LIBFLAME_PATH"
    )
endif()

# Assign the include directory based on the resolved root
set(LIBFLAME_INCLUDE_DIR "${DETERMINED_ROOT}/include" CACHE PATH "Path to libflame headers")

# Search for the static or shared library file within the root's lib/ directory
find_library(LIBFLAME_LIBRARY
    NAMES flame libflame
    PATHS "${DETERMINED_ROOT}/lib"
          "${DETERMINED_ROOT}/lib64"
    NO_DEFAULT_PATH
    REQUIRED
)

message(STATUS "Using LIBFLAME_INCLUDE_DIR: ${LIBFLAME_INCLUDE_DIR}")
message(STATUS "Using LIBFLAME_LIBRARY: ${LIBFLAME_LIBRARY}")

# ---------------------------------------------------------------------------
# Build test executable
# ---------------------------------------------------------------------------
file(GLOB SRC_FILES "src/*.c")
add_executable(test_libflame.x  ${SRC_FILES})

# Include headers for both libflame and the local test directory definitions
target_include_directories(test_libflame.x PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${LIBFLAME_INCLUDE_DIR}
)

# Link everything together (Math library, libflame, and BLAS targets)
target_link_libraries(test_libflame.x PRIVATE
    ${LIBFLAME_LIBRARY}
    ${BLAS_LIBRARIES}
    m
)

# ---------------------------------------------------------------------------
# Add test command
# ---------------------------------------------------------------------------
# Copy the operational configuration files to the build directory for execution
configure_file(input.global.general input.global.general COPYONLY)
configure_file(input.global.operations input.global.operations COPYONLY)

enable_testing()

# Define the test case
add_test(NAME test_libflame
         COMMAND test_libflame.x
)

set_tests_properties(test_libflame PROPERTIES
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)
