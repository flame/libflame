/******************************************************************************
* Copyright (C) 2021, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file libflame_test.hh
 *  libflame_test.hh defines all the common functions which is required
 *  to validate CPP template interfaces for all libflame modules
 *  */

#ifndef LIBFLAME_TEST_HH
#define LIBFLAME_TEST_HH

#include "libflame_interface.hh"

#ifndef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
#define FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
#endif

/*! @brief computeDiff is function template to compare contents of 2 buffers.
			T can be float, double.
 * @details
 * \b Purpose:
    \verbatim
	  computeDiff is function template to compare contents of 2 buffers.
    T can be float, double.
    \endverbatim
 * @param[in] size
          size is integer.
          Size of the input buffers.
 * @param[in] Out
          Out is REAL or DOUBLE PRECISION array.
          Contents of Buffer.
 * @param[in] Out_ref
          Out_ref is REAL or DOUBLE PRECISION array.
          Contents of Reference buffer.
 * @return DOUBLE
          Returns differences value after comparing output of C and CPP based
          library APIs.
 * */
template<typename T>
double computeDiff(int size, T *Out, T *Out_ref)
{
  int index;
  T diff = 0;

  for (index = 0; index < size; index ++) {
    diff += abs(Out_ref[index] - Out[index]);
  }
  return diff;
}

/*! @brief computeErrorComplex is function template to compare contents of
      2 buffers. T can be scomplex, dcomplex.
 * @details
 * \b Purpose:
    \verbatim
	  computeErrorComplex is function template to compare contents of 2 buffers.
    T can be scomplex, dcomplex.
    \endverbatim
 * @param[in] size
          size is integer.
          Size of the input buffers.
 * @param[in] Out
          Out is COMPLEX or COMPLEX*16 array.
          Contents of Buffer.
 * @param[in] Out_ref
          Out_ref is COMPLEX or COMPLEX*16 array.
          Contents of Reference buffer.
 * @return DOUBLE
          Returns differences value after comparing output of C and CPP based
          library APIs.
 * */
template<typename T>
double computeErrorComplex(int size, T *Out, T *Out_ref)
{
  int index;
  double diff = 0;
  for (index = 0; index < size; index ++) {
    diff += abs(Out_ref[index].real - Out[index].real);
    diff += abs(Out_ref[index].imag - Out[index].imag);
  }
  return diff;
}

/*! @brief computeError is function template to compare contents of 2 buffers.
			T can be float, double, scomplex, dcomplex.
 * @details
 * \b Purpose:
    \verbatim
	  computeError is function template to compare contents of 2 buffers.
    T can be float, double, scomplex, dcomplex.
    \endverbatim
 * @param[in] n
          n is integer.
          n of the number of rows.
 * @param[in] m
          m is integer.
          m of the number of columns.
          m*n gives size of the array/buffer.
 * @param[in] Out
          Out is REAL or DOUBLE PRECISION or COMPLEX or COMPLEX*16 array.
          Contents of Buffer.
 * @param[in] Out_ref
          Out_ref is REAL or DOUBLE PRECISION or COMPLEX or COMPLEX*16 array.
          Contents of Reference buffer.
 * @return DOUBLE
          Returns differences value after comparing output of C and CPP based
          library APIs.
 * */
template<typename T>
double computeError(int n, int m, T *out, T *out_ref)
{
  double retDiff = 0;
  if (typeid(T) == typeid(float)) {
    retDiff = computeDiff<float>(n*m, (float*)out, (float*)out_ref);
  } else if (typeid(T) == typeid(double)) {
    retDiff = computeDiff<double>(n*m, (double *)out, (double *)out_ref);
  } else if (typeid(T) == typeid(int)) {
    retDiff = computeDiff<int>(n*m, (int *)out, (int *)out_ref);
  } else if (typeid(T) == typeid(scomplex)) {
    retDiff = computeErrorComplex<scomplex>(n*m, (scomplex *)out,
                                            (scomplex *)out_ref);
  } else if (typeid(T) == typeid(dcomplex)) {
    retDiff = computeErrorComplex<dcomplex>(n*m, (dcomplex *)out,
                                            (dcomplex *)out_ref);
  }
  return retDiff;
}

/*! @brief allocate_init_buffer is function to allocate memory for buffer
      and initialise memory with random values.

 * @details
 * \b Purpose:
    \verbatim
    allocate_init_buffer is function to allocate memory for buffer
    and initialise memory with random values.
    \endverbatim
 * @param[in] aIn
          aIn is integer array.
          Buffer pointer to allocate and initialize based on size.
 * @param[in] aRef
          aRef is integer array.
          Reference buffer pointer to allocate and initialize based on size.
 * @param[in] size
          size is integer.
          Size of the input buffers.
 * @return VOID
          Nothing.
 * */
void allocate_init_buffer(int *&aIn, int *&aRef, int size)
{
  aIn = new int [size];
  aRef = new int [size];
  for (int index = 0; index < size; index++) {
    aIn[index] = ((int) rand() / ((int) RAND_MAX / 2.0)) - 1.0;
    aRef[index] =  aIn[index];
  }
}

/*! @brief allocate_init_buffer is function to allocate memory for buffer
      and initialise memory with random values.

 * @details
 * \b Purpose:
    \verbatim
    allocate_init_buffer is function to allocate memory for buffer
    and initialise memory with random values.
    \endverbatim
 * @param[in] aIn
          aIn is REAL array.
          Buffer pointer to allocate and initialize based on size.
 * @param[in] aRef
          aRef is REAL array.
          Reference buffer pointer to allocate and initialize based on size.
 * @param[in] size
          size is integer.
          Size of the input buffers.
 * @return VOID
          Nothing.
 * */
void allocate_init_buffer(float *&aIn, float *&aRef, int size)
{
  aIn = new float [size];
  aRef = new float [size];
  for (int index = 0; index < size; index++) {
    aIn[index] = ((float) rand() / ((float) RAND_MAX / 2.0));
    aRef[index] =  aIn[index] ;
  }
}

/*! @brief allocate_init_buffer is function to allocate memory for buffer
      and initialise memory with random values.

 * @details
 * \b Purpose:
    \verbatim
    allocate_init_buffer is function to allocate memory for buffer
    and initialise memory with random values.
    \endverbatim
 * @param[in] aIn
          aIn is DOUBLE PRECISION array.
          Buffer pointer to allocate and initialize based on size.
 * @param[in] aRef
          aRef is DOUBLE PRECISION array.
          Reference buffer pointer to allocate and initialize based on size.
 * @param[in] size
          size is integer.
          Size of the input buffers.
 * @return VOID
          Nothing.
 * */
void allocate_init_buffer(double *&aIn, double *&aRef, int size)
{
  aIn = new double [size];
  aRef = new double [size];
  for (int index = 0; index < size; index++) {
    aIn[index] = ((double) rand() / ((double) RAND_MAX / 2.0));
    aRef[index] =  aIn[index] ;
  }
}

/*! @brief allocate_init_buffer is function to allocate memory for buffer
      and initialise memory with random values.

 * @details
 * \b Purpose:
    \verbatim
    allocate_init_buffer is function to allocate memory for buffer
    and initialise memory with random values.
    \endverbatim
 * @param[in] aIn
          aIn is COMPLEX array.
          Buffer pointer to allocate and initialize based on size.
 * @param[in] aRef
          aRef is COMPLEX array.
          Reference buffer pointer to allocate and initialize based on size.
 * @param[in] size
          size is integer.
          Size of the input buffers.
 * @return VOID
          Nothing.
 * */
void allocate_init_buffer(scomplex *&aIn, scomplex *&aRef, int size)
{
  aIn = new scomplex [size];
  aRef = new scomplex [size];
  for (int index = 0; index < size; index++) {
    float real = ((float) rand() / ((float) RAND_MAX / 2.0)) - 1.0;
    float imag = ((float) rand() / ((float) RAND_MAX / 2.0)) - 1.0;
    aIn[index].real = real;
    aIn[index].imag = imag;
    aRef[index] = aIn[index];
  }
}

/*! @brief allocate_init_buffer is function to allocate memory for buffer
      and initialise memory with random values.

 * @details
 * \b Purpose:
    \verbatim
    allocate_init_buffer is function to allocate memory for buffer
    and initialise memory with random values.
    \endverbatim
 * @param[in] aIn
          aIn is COMPLEX*16 array.
          Buffer pointer to allocate and initialize based on size.
 * @param[in] aRef
          aRef is COMPLEX*16 array.
          Reference buffer pointer to allocate and initialize based on size.
 * @param[in] size
          size is integer.
          Size of the input buffers.
 * @return VOID
          Nothing.
 * */
void allocate_init_buffer(dcomplex *&aIn, dcomplex *&aRef, int size)
{
  aIn = new dcomplex [size];
  aRef = new dcomplex [size];

  for (int index = 0; index < size; index++) {
    double real = ((double) rand() / ((double) RAND_MAX / 2.0)) - 1.0;
    double imag = ((double) rand() / ((double) RAND_MAX / 2.0)) - 1.0;
    aIn[index].real = real;
    aIn[index].imag = imag;
    aRef[index] = aIn[index];
  }
}

/*! @brief allocate_init_buffer is function to allocate memory for buffer
      and initialise memory with specified values.

 * @details
 * \b Purpose:
    \verbatim
    allocate_init_buffer is function to allocate memory for buffer
    and initialise memory with specified values.
    \endverbatim
 * @param[in] aIn
          aIn is integer array.
          Buffer pointer to allocate and initialize based on size.
 * @param[in] aRef
          aRef is integer array.
          Reference buffer pointer to allocate and initialize based on size.
 * @param[in] size
          size is integer.
          Size of the input buffers.
 * @param[in] value
          value is integer.
          value to initialize in the input buffers.
 * @return VOID
          Nothing.
 * */
void allocate_init_buffer(int *&aIn, int *&aRef, int size, int value)
{
  aIn =  new int [size];
  aRef = new int [size];
  for (int index = 0; index < size; index++) {
    aIn[index] = value;
    aRef[index] = value;
  }
}

/*! @brief allocate_init_buffer is function to allocate memory for buffer
      and initialise memory with specified values.

 * @details
 * \b Purpose:
    \verbatim
    allocate_init_buffer is function to allocate memory for buffer
    and initialise memory with specified values.
    \endverbatim
 * @param[in] aIn
          aIn is REAL array.
          Buffer pointer to allocate and initialize based on size.
 * @param[in] aRef
          aRef is REAL array.
          Reference buffer pointer to allocate and initialize based on size.
 * @param[in] size
          size is integer.
          Size of the input buffers.
 * @param[in] value
          value is integer.
          value to initialize in the input buffers.
 * @return VOID
          Nothing.
 * */
void allocate_init_buffer(float *&aIn, float *&aRef, int size, int value)
{
  aIn = new float [size];
  aRef = new float [size];
  for (int index = 0; index < size; index++) {
    aIn[index] = value;
    aRef[index] = value;
  }
}

/*! @brief allocate_init_buffer is function to allocate memory for buffer
      and initialise memory with specified values.

 * @details
 * \b Purpose:
    \verbatim
    allocate_init_buffer is function to allocate memory for buffer
    and initialise memory with specified values.
    \endverbatim
 * @param[in] aIn
          aIn is DOUBLE PRECISION array.
          Buffer pointer to allocate and initialize based on size.
 * @param[in] aRef
          aRef is DOUBLE PRECISION array.
          Reference buffer pointer to allocate and initialize based on size.
 * @param[in] size
          size is integer.
          Size of the input buffers.
 * @param[in] value
          value is integer.
          value to initialize in the input buffers.
 * @return VOID
          Nothing.
 * */
void allocate_init_buffer(double *&aIn, double *&aRef, int size, int value)
{
  aIn = new double [size];
  aRef = new double [size];
  for (int index = 0; index < size; index++) {
    aIn[index] = value;
    aRef[index] = value;
  }
}

/*! @brief allocate_init_buffer is function to allocate memory for buffer
      and initialise memory with specified values.

 * @details
 * \b Purpose:
    \verbatim
    allocate_init_buffer is function to allocate memory for buffer
    and initialise memory with specified values.
    \endverbatim
 * @param[in] aIn
          aIn is COMPLEX array.
          Buffer pointer to allocate and initialize based on size.
 * @param[in] aRef
          aRef is COMPLEX array.
          Reference buffer pointer to allocate and initialize based on size.
 * @param[in] size
          size is integer.
          Size of the input buffers.
 * @param[in] value
          value is integer.
          value to initialize in the input buffers.
 * @return VOID
          Nothing.
 * */
void allocate_init_buffer(scomplex *&aIn, scomplex *&aRef, int size, int value)
{
  aIn = new scomplex [size];
  aRef = new scomplex [size];
  for (int index = 0; index < size; index++) {
    aIn[index].imag = value;
    aIn[index].real = value;
    aRef[index] = aIn[index];
  }
}

/*! @brief allocate_init_buffer is function to allocate memory for buffer
      and initialise memory with specified values.

 * @details
 * \b Purpose:
    \verbatim
    allocate_init_buffer is function to allocate memory for buffer
    and initialise memory with specified values.
    \endverbatim
 * @param[in] aIn
          aIn is COMPLEX*16 array.
          Buffer pointer to allocate and initialize based on size.
 * @param[in] aRef
          aRef is COMPLEX*16 array.
          Reference buffer pointer to allocate and initialize based on size.
 * @param[in] size
          size is integer.
          Size of the input buffers.
 * @param[in] value
          value is integer.
          value to initialize in the input buffers.
 * @return VOID
          Nothing.
 * */
void allocate_init_buffer(dcomplex *&aIn, dcomplex *&aRef, int size, int value)
{
  aIn = new dcomplex [size];
  aRef = new dcomplex [size];
  for (int index = 0; index < size; index++) {
    aIn[index].imag = value;
    aIn[index].real = value;
    aRef[index] = aIn[index];
  }
}

/*! @brief print_array_data is function to print the contents of array
      If MAX_ARRAY_PRINT_SIZE is defined, then all the contents of
      array will be printed. Else if ARRAY_PRINT_SIZE is defined and 
      it is less than the arraysize parameter, then the contents with
      ARRAY_PRINT_SIZE will be printed. Otherwise,
      the contents of array will be printed till arraysize.

 * @details
 * \b Purpose:
    \verbatim
    print_array_data is function to print the contents of array
    If MAX_ARRAY_PRINT_SIZE is defined, then all the contents of
    array will be printed. Else if ARRAY_PRINT_SIZE is defined and 
    it is less than the arraysize parameter, then the contents with
    ARRAY_PRINT_SIZE will be printed. Otherwise,
    the contents of array will be printed till arraysize.
    \endverbatim
 * @param[in] arrayname
          arrayname is CHAR array.
          Contains the name of the array.
 * @param[in] aray
          array is REAL array.
          Array pointer to print the contents.
 * @param[in] arraysize
          arraysize is integer.
          Size of the input array.

 * @return VOID
           Nothing.
 * */
template <typename T>
void print_array_data(char *arrayname, T *array, integer arraysize)
{
  integer size = arraysize;

  #if defined(MAX_ARRAY_PRINT_SIZE) && (MAX_ARRAY_PRINT_SIZE == 0) && \
      defined(ARRAY_PRINT_SIZE) && (ARRAY_PRINT_SIZE > 0)
      if (ARRAY_PRINT_SIZE < arraysize) {
        size = ARRAY_PRINT_SIZE;
      }
  #endif
  // Specifier default is to print float, changes based on integer.
  char specifier[] = "%f ";
  if (typeid(T) == typeid(integer)) {
    specifier[1] = 'd';
  }
  PRINTF("Printing %s array with %d size.\n", arrayname, size);
  for (int index = 0; index < size; index++) {
    PRINTF(specifier, array[index]);
  }
  PRINTF("\n");
}

/*! @brief print_array_complex_data is function to print the contents of array
      If MAX_ARRAY_PRINT_SIZE is defined, then all the contents of
      array will be printed. Else if ARRAY_PRINT_SIZE is defined and 
      it is less than the arraysize parameter, then the contents with
      ARRAY_PRINT_SIZE will be printed. Otherwise,
      the contents of array will be printed till arraysize.

 * @details
 * \b Purpose:
    \verbatim
    print_array_complex_data is function to print the contents of array
    If MAX_ARRAY_PRINT_SIZE is defined, then all the contents of
    array will be printed. Else if ARRAY_PRINT_SIZE is defined and 
    it is less than the arraysize parameter, then the contents with
    ARRAY_PRINT_SIZE will be printed. Otherwise,
    the contents of array will be printed till arraysize.
    \endverbatim
 * @param[in] arrayname
          arrayname is CHAR array.
          Contains the name of the array.
 * @param[in] array
          array is REAL array.
          Array pointer to print the contents.
 * @param[in] arraysize
          arraysize is integer.
          Size of the input array.

 * @return VOID
           Nothing.
 * */
template <typename T>
void print_array_complex_data(char *arrayname, T *array, integer arraysize)
{
  integer size = arraysize;
  
  #if defined(MAX_ARRAY_PRINT_SIZE) && (MAX_ARRAY_PRINT_SIZE == 0) && \
      defined(ARRAY_PRINT_SIZE) && (ARRAY_PRINT_SIZE > 0)
      if (ARRAY_PRINT_SIZE < arraysize) {
        size = ARRAY_PRINT_SIZE;
      }
  #endif
  PRINTF("Printing %s array with %d size.\n", arrayname, size);
  for (int index = 0; index < size; index++) {
    PRINTF("%f+%fi ", array[index].real, array[index].imag);
  }
  PRINTF("\n");
}

/*! @brief print_array is function to print the contents of array.
      Based on the template data typename, it will call respective function.
      print_array_data() will be called to print the contents of float, 
      double and integer arrays.
      print_array_complex_data() will be called to print the contents of
      scomplex, dcomplex arrays.

 * @details
 * \b Purpose:
    \verbatim
    print_array is function to print the contents of array.
    Based on the template data typename, it will call respective function.
    print_array_data() will be called to print the contents of float, 
    double and integer arrays.
    print_array_complex_data() will be called to print the contents of
    scomplex, dcomplex arrays.
    \endverbatim
 * @param[in] arrayname
          arrayname is CHAR array.
          Contains the name of the array.
 * @param[in] array
          array is REAL array.
          Array pointer to print the contents.
 * @param[in] arraysize
          arraysize is integer.
          Size of the input array.

 * @return VOID
           Nothing.
 * */
template <typename T>
void print_array(char *arrayname, T *array, integer arraysize)
{
  if (typeid(T) == typeid(float)) {
    print_array_data<float>(arrayname, (float *)array, arraysize);
  } else if (typeid(T) == typeid(double)) {
    print_array_data<double>(arrayname, (double *)array, arraysize);
  } else if (typeid(T) == typeid(integer)) {
    print_array_data<integer>(arrayname, (integer *)array, arraysize);
  } else if (typeid(T) == typeid(scomplex)) {
    print_array_complex_data<scomplex>(arrayname, (scomplex *)array, arraysize);
  } else if (typeid(T) == typeid(dcomplex)) {
    print_array_complex_data<dcomplex>(arrayname, (dcomplex *)array, arraysize);
  }
}

#endif        //  #ifndef LIBFLAME_TEST_HH