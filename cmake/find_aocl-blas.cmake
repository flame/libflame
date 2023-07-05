# ########################################################################
#Copyright(c) 2023 Advanced Micro Devices, Inc.
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files(the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and / or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following                       conditions:
#
#The above copyright notice and this permission notice shall be included in
#all copies or substantial portions of the                               Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#THE SOFTWARE.
#
# ########################################################################

# ============= aocl function ================
function(aocl_libs)

  IF(FLA_ENABLE_ILP64) 
    SET(ILP_DIR "ILP64")
  ELSE(FLA_ENABLE_ILP64)
    SET(ILP_DIR "LP64")
  ENDIF(FLA_ENABLE_ILP64)

  IF(WIN32)
    SET(CMAKE_FIND_LIBRARY_PREFIXES "")
    SET(CMAKE_FIND_LIBRARY_SUFFIXES ".lib") 
    find_library(AOCL_BLAS_LIB
    NAMES AOCL-LibBlis-Win-MT AOCL-LibBlis-Win-MT-dll AOCL-LibBlis-Win AOCL-LibBlis-Win-dll
    HINTS ${AOCL_ROOT}/blis ${AOCL_ROOT}/amd-blis ${AOCL_ROOT}
    PATH_SUFFIXES "lib/${ILP_DIR}" "lib_${ILP_DIR}" "lib"
    DOC "AOCL-BLAS library"
    )

    #====Headers
    find_path(AOCL_BLAS_INCLUDE_DIR
    NAMES blis.h  cblas.h
    HINTS ${AOCL_ROOT}/amd-blis ${AOCL_ROOT}/blis ${AOCL_ROOT}
    PATH_SUFFIXES "include/${ILP_DIR}" "include_${ILP_DIR}" "include" "include/blis"
    DOC "AOCL-BLAS headers"
    )

  ELSE(WIN32)   
    SET(CMAKE_FIND_LIBRARY_PREFIXES "lib")
    IF(BUILD_SHARED_LIBS)
      SET(CMAKE_FIND_LIBRARY_SUFFIXES ".so")
    ELSE(BUILD_SHARED_LIBS)
      SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    ENDIF(BUILD_SHARED_LIBS) 

    find_library(AOCL_BLAS_LIB
    NAMES blis-mt blis
    HINTS ${AOCL_ROOT}/blis ${AOCL_ROOT}/amd-blis ${AOCL_ROOT}
    PATH_SUFFIXES "lib/${ILP_DIR}" "lib_${ILP_DIR}" "lib"
    DOC "AOCL-BLAS library"
    )
    
    #====Headers
    find_path(AOCL_BLAS_INCLUDE_DIR
    NAMES blis.h  blis.hh  cblas.h  cblas.hh
    HINTS ${AOCL_ROOT}/blis ${AOCL_ROOT}/amd-blis ${AOCL_ROOT}
    PATH_SUFFIXES "include/${ILP_DIR}" "include_${ILP_DIR}" "include" "include/blis"
    DOC "AOCL-BLAS headers"
    )

  ENDIF(WIN32)

  #===========
  if(AOCL_BLAS_LIB AND AOCL_BLAS_INCLUDE_DIR)
    set(AOCL_BLAS_FOUND true PARENT_SCOPE)
    set(BLAS_LIBRARY ${AOCL_BLAS_LIB})  
  else()    
    message (FATAL_ERROR "Error: could not find a suitable installation of AOCL-BLAS Library in \$AOCL_ROOT=${AOCL_ROOT}")
  endif()

  set(BLAS_LIBRARY ${BLAS_LIBRARY} PARENT_SCOPE)

endfunction(aocl_libs)

#==================main=================
# clear to avoid endless appending on subsequent calls
set(BLAS_LIBRARY)
unset(BLAS_INCLUDE_DIR)

if(DEFINED ENV{AOCL_ROOT})            
    SET(AOCL_ROOT $ENV{AOCL_ROOT})
    message(STATUS "AOCL_ROOT set via environment variable is ${AOCL_ROOT}")
    if(NOT EXISTS ${AOCL_ROOT})
			message(FATAL_ERROR "\n Invalid path to AOCL_ROOT \n")
		endif()
elseif(AOCL_ROOT)
      SET(AOCL_ROOT ${AOCL_ROOT})
      message(STATUS "AOCL_ROOT set from cmake option is ${AOCL_ROOT}")
      if(NOT EXISTS ${AOCL_ROOT})
        message(FATAL_ERROR "\n Invalid path to AOCL_ROOT \n")
      endif()
else()      
      set(AOCL_ROOT ${CMAKE_INSTALL_PREFIX})
endif()

aocl_libs()

set(BLAS_LIBRARY ${BLAS_LIBRARY})
set(BLAS_INCLUDE_DIR ${AOCL_BLAS_INCLUDE_DIR})

message(STATUS "AOCL-BLAS LIBRARY= ${BLAS_LIBRARY}")
message(STATUS "AOCL-BLAS INCLUDE DIRS= ${BLAS_INCLUDE_DIR}")

mark_as_advanced(BLAS_LIBRARY BLAS_INCLUDE_DIR)
