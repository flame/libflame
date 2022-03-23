/******************************************************************************
* Copyright (C) 2021-2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file libflame_test.hh
 *  libflame_test.hh defines all the common functions which is required
 *  to validate CPP template interfaces for all libflame modules
 *  */

#ifndef LIBFLAME_TEST_HH
#define LIBFLAME_TEST_HH

#include "libflame_interface.hh"
#include "libflame_utils_interface.hh"
#include <sys/time.h>
#include <time.h>

#ifndef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
#define FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
#endif

// Global variables.
double ref_time_sec = 0.0;
float s_zero = 0, s_one = 1, s_n_one = -1;
double d_zero = 0, d_one = 1, d_n_one = -1;
scomplex c_zero = {0,0}, c_one = {1,0}, c_n_one = {-1,0};
dcomplex z_zero = {0,0}, z_one = {1,0}, z_n_one = {-1,0};
void *n_one = NULL, *one = NULL, *zero = NULL;

#define DRAND()  ( ( double ) rand() / ( ( double ) RAND_MAX / 2.0F ) ) - 1.0F;
#define SRAND()  ( float ) ( ( double ) rand() / ( ( double ) RAND_MAX / 2.0F ) ) - 1.0F;

template <typename T>
void alpha_beta_init()
{
  if (typeid(T) == typeid(float)) {
    n_one = &s_n_one;
    one = &s_one;
    zero = &s_zero;
  } else if (typeid(T) == typeid(double)) {
    n_one = &d_n_one;
    one = &d_one;
    zero = &d_zero;
  } else if (typeid(T) == typeid(scomplex)) {
    n_one = &c_n_one;
    one = &c_one;
    zero = &c_zero;
  } else if (typeid(T) == typeid(dcomplex)) {
    n_one = &z_n_one;
    one = &z_one;
    zero = &z_zero;
  }
}

/* create vector of given datatype*/
template <typename T>
void create_vector(T **A, integer M)
{
  *A = (T *)malloc(M * sizeof(T));

  if (*A == NULL) {
    fprintf( stderr, "malloc() returned NULL pointer\n");
    abort();
  }
}

template <typename T>
void create_realtype_vector(T **A, integer M)
{
  *A = NULL;

  if((typeid(T) == typeid(float)) || (typeid(T) == typeid(scomplex))) {
    *A = (float *)malloc(M * sizeof(float));
  } else {
    *A = (double *)malloc(M * sizeof(double));
  }
  if (*A == NULL) {
    fprintf( stderr, "malloc() returned NULL pointer\n");
    abort();
  }
}

/* free vector */
void free_vector(void *A)
{
  if(A) {
    free(A);
  }
}

/* Initialize vector with random values */
template <typename T>
void rand_vector(T *A, integer M, integer LDA)
{
  integer i;

  if (typeid(T) == typeid(float)) {
    for (i = 0; i < M; i++) {
      ((float *)A)[i * LDA] = SRAND();
    }
  } else if (typeid(T) == typeid(double)) {
    for (i = 0; i < M; i++) {
      ((double *)A)[i * LDA] = SRAND();
    }
  } else if (typeid(T) == typeid(scomplex)) {
    for (i = 0; i < M; i++) {
      ((scomplex *)A)[i * LDA].real = SRAND();
      ((scomplex *)A)[i * LDA].imag = SRAND();
    }
  } else if (typeid(T) == typeid(dcomplex)) {
    for (i = 0; i < M; i++) {
      ((dcomplex *)A)[i * LDA].real = SRAND();
      ((dcomplex *)A)[i * LDA].imag = SRAND();
    }
  }
}

/* Copy a vector */
template <typename T>
void copy_vector(integer M, T *A, integer LDA, T *B, integer LDB)
{
  if (typeid(T) == typeid(integer)) {
    integer i;

    for (i = 0; i < M; i++) {
      ((T *)B)[ i * LDB ] = ((T *)A)[ i * LDA ];
    }
  } else {
    libflame_utils::copy<T>(&M, A, &LDA, B, &LDB);
  }
}

template <typename T>
void copy_realtype_vector(integer M, T *A, integer LDA, T *B, integer LDB)
{
  if((typeid(T) == typeid(float)) || (typeid(T) == typeid(scomplex))) {
    libflame_utils::copy<float>(&M, A, &LDA, B, &LDB);
  } else {
    libflame_utils::copy<double>(&M, A, &LDA, B, &LDB);
  }
}

/* create matrix of given datatype*/
template <typename T>
void create_matrix(T **A, integer M, integer N)
{
  *A = (T *)malloc(M * N * sizeof(T));
  if (*A == NULL) {
    fprintf (stderr, "malloc() returned NULL pointer\n");
    abort();
  }
}

template <typename T>
void create_realtype_matrix(T **A, integer M, integer N)
{
  *A = NULL;

  if ((typeid(T) == typeid(float)) || (typeid(T) == typeid(scomplex))) {
    *A = (float *)malloc(M * N * sizeof(float));
  } else {
    *A = (double *)malloc(M * N * sizeof(double));
  }

  if (*A == NULL) {
    fprintf (stderr, "malloc() returned NULL pointer\n");
    abort();
  }
}

template <typename T>
T* get_m_ptr(T *A, integer M, integer N, integer LDA)
{
  T *mat = ((T *)A) + M + N * LDA;
  return mat;
}

/* free matrix */
void free_matrix(void *A)
{
  if (A) {
    free(A);
  }
}

template <typename T>
void rand_matrix_sd(T *A, integer M, integer N, integer LDA)
{
  integer i, j;
  for (i = 0; i < N; i++) {
    for (j = 0; j < M; j++) {
      ((T *)A)[i * LDA + j] = SRAND();
    }
  }
}

template <typename T>
void rand_matrix_cz(T *A, integer M, integer N, integer LDA)
{
  integer i, j;
  for(i = 0; i < N; i++) {
    for(j = 0; j < M; j++) {
      ((T *)A)[i * LDA + j].real = SRAND();
      ((T *)A)[i * LDA + j].imag = SRAND();
    }
  }
}

/* Initialize matrix with random values */
template <typename T>
void rand_matrix(T *A, integer M, integer N, integer LDA)
{
  if (typeid(T) == typeid(float)) {
    rand_matrix_sd<float>((float *)A, M, N, LDA);
  } else if (typeid(T) == typeid(double)) {
    rand_matrix_sd<double>((double *)A, M, N, LDA);
  } else if (typeid(T) == typeid(scomplex)) {
    rand_matrix_cz<scomplex>((scomplex *)A, M, N, LDA);
  } else if (typeid(T) == typeid(dcomplex)) {
    rand_matrix_cz<dcomplex>((dcomplex *)A, M, N, LDA);
  }
}

/* Initialize symmetric matrix with random values */
template <typename T>
void rand_sym_matrix(T *A, integer M, integer N, integer LDA)
{
  integer i, j;

  if (typeid(T) == typeid(float)) {
    for (i = 0; i < N; i++ ) {
      for (j = i; j < M; j++) {
        ((T *)A)[i * LDA + j] = SRAND();
        ((T *)A)[j * LDA + i] = ((T *)A)[i * LDA + j];
      }
    }
  } else if (typeid(T) == typeid(double)) {
    for (i = 0; i < N; i++) {
      for (j = i; j < M; j++) {
        ((T *)A)[i * LDA + j] = SRAND();
        ((T *)A)[j * LDA + i] = ((T *)A)[i * LDA + j];
      }
    }
  }
  else if ((typeid(T) == typeid(scomplex)) || 
          (typeid(T) == typeid(dcomplex))) {
    for (i = 0; i < N; i++) {
      for (j = i; j < M; j++) {
        ((T *)A)[i * LDA + j].real = SRAND();
        ((T *)A)[i * LDA + j].imag = SRAND();
      }
    }

    for (i = 0; i < N; i++) {
      for (j = i; j < M; j++) {
        ((T *)A)[j * LDA + i].real = ((T *)A)[i * LDA + j].real;
        ((T *)A)[j * LDA + i].imag = ((T *)A)[i * LDA + j].imag;
      }
    }
  }
}

/* Copy a matrix */
template <typename T>
void copy_matrix(char *uplo, integer M, integer N, T *A, integer LDA, 
                 T *B, integer LDB)
{
  if (typeid(T) == typeid(integer)) {
    integer i, j;

    for (i = 0; i < N; i++ ) {
      for (j = 0; j < M; j++) {
        ((integer *)B)[ i * LDB + j ] = ((integer *)A)[ i * LDA + j ];
      }
    }
  }
  else {
    libflame::lacpy<T>(uplo, &M, &N, A, &LDA, B, &LDB);
  }
}

template <typename T>
void copy_realtype_matrix(char *uplo, integer M, integer N, T *A, integer LDA,
                          T *B, integer LDB)
{
  if ((typeid(T) == typeid(float)) || (typeid(T) == typeid(scomplex))) {
    slacpy_(uplo, &M, &N, A, &LDA, B, &LDB);
  } else {
    dlacpy_(uplo, &M, &N, A, &LDA, B, &LDB);
  }
}

/* Initialize a matrix with zeros */
template <typename T>
void reset_matrix(integer M, integer N, T *A, integer LDA)
{
  integer i, j;

  if (typeid(T) == typeid(integer)) {
    for (i = 0; i < N; i++ ) {
      for (j = 0; j < M; j++) {
        ((integer *)A)[ i * LDA + j ] = 0.f;
      }
    }
  } else {
    alpha_beta_init<T>();
    libflame::laset<T>("A", &M, &N, (T *)zero, (T *)zero, A, &LDA);
  }
}

/* Set a matrix to identity */
template <typename T>
void set_identity_matrix(integer M, integer N, T *A, integer LDA)
{
  alpha_beta_init<T>();
  libflame::laset<T>("A", &M, &N, (T *)zero, (T *)one, A, &LDA);
}

void z_div_t(dcomplex *cp, dcomplex *ap, dcomplex *bp)
{
  dcomplex a = *ap;
  dcomplex b = *bp;
  double temp;

  temp = b.real * b.real + b.imag * b.imag;
  if (!temp) {
    fprintf (stderr, "z_div_t : temp is zero. Abort\n");
    abort();
  }

  cp->real = (a.real * b.real + a.imag * b.imag) / temp;
  cp->imag = (a.imag * b.real - a.real * b.imag) / temp;
}


/* Division of complex types */
void c_div_t(scomplex *cp, scomplex *ap, scomplex *bp)
{
  scomplex a = *ap;
  scomplex b = *bp;
  float temp;

  temp = b.real * b.real + b.imag * b.imag;
  if (!temp) {
    fprintf (stderr, "z_div_t : temp is zero. Abort\n");
    abort();
  }

  cp->real = (a.real * b.real + a.imag * b.imag) / temp;
  cp->imag = (a.imag * b.real - a.real * b.imag) / temp;
}

/* work value calculation */
template <typename T>
integer get_work_value(T *work)
{
  if (!work) {
    return 0;
  }

  if ((typeid(T) == typeid(float)) || (typeid(T) == typeid(scomplex))) {
    return (integer) (*(float*)work);
  } else {
    return (integer) (*(double*)work);
  }
}

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

/*! @brief  getrf_internal() function template calls C and CPP based 
      GETRF APIs with valid test values and returns the outout buffers
      and info.
      T can be float, double, scomplex, dcomplex.
 * @details
 * \b Purpose:
    \verbatim
    getrf_internal is function template for getrf() functions.
    T can be float, double, scomplex, dcomplex.
    
    getrf_internal() function template calls C and CPP based GETRF APIs with
    valid test values and returns the outout buffers and info.
    abuff, ipivbuff, arefbuff, ipivrefbuff, info_cpp, info_ref will be
    updated as output.
    if abuff, ipivbuff buffers are NULL, then respective CPP
    function call is not performed.
    if arefbuff, ipivrefbuff buffers are NULL, then respective C
    function call is not performed.
    NOTE: It is caller's responsibility to free the buffers after the function
    call.
    \endverbatim
  
 * @params[in] m
          M is INTEGER
          The number of rows of the matrix A.  M >= 0.
 * @params[in] n
          N is INTEGER
          The number of columns of the matrix A.  N >= 0.
 * @params[in] abuff
          A is REAL array, dimension (LDA,N)
          On entry, the M-by-N matrix to be factored.
          On exit, the factors L and U from the factorization
          A = P*L*U; the unit diagonal elements of L are not stored.
 * @params[in] lda
          LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,M).
 * @params[in] ipivbuff
          IPIV is INTEGER array, dimension (min(M,N))
          The pivot indices; for 1 <= i <= min(M,N), row i of the
          matrix was interchanged with row IPIV(i).
 * @params[in] arefbuff
          arefbuff is REAL array, dimension (LDA,N)
          arefbuff is similar to abuff & used for reference C function.
 * @params[in] ipivrefbuff
          ipivrefbuff is INTEGER array, dimension (min(M,N))
          ipivrefbuff is similar to ipivbuff & used for reference C function.
 * @params[out] info_cpp
          info_cpp is integer to get INFO returned by CPP function.
 * @params[out] info_ref
          info_ref is integer to get INFO returned by reference C function.

 * @return VOID
           Nothing.
 * */
template<typename T>
void getrf_internal(integer m, integer n, T* abuff,
        integer lda, integer* ipivbuff, T* arefbuff,
        integer* ipivrefbuff, integer* info_cpp,
        integer* info_ref)
{
  typedef int (*fptr_NL_LAPACK_getrf)(integer* m, integer* n, T* a,
                  integer* lda, integer* ipiv, integer* info);
  fptr_NL_LAPACK_getrf getrf_ref = NULL;
  
  PRINTF("%s() Entry...\n", __FUNCTION__);
  if ((m <= 0) || \
      (n <= 0) || \
      (lda <= 0) || \
      (info_cpp == NULL) || \
      (info_ref == NULL)) {
    PRINTF("Please pass the valid inputs/buffers. Returning from function.\n");
    return;
  }
  
  if (((abuff == NULL) && (arefbuff == NULL)) ||
      ((ipivbuff == NULL) && (ipivrefbuff == NULL))) {
    PRINTF("Please pass the valid buffers. Returning from function.\n");
    return;
  }
  
  // Print input values other than arrays.
  #if (defined(PRINT_INPUT_VALUES) && (PRINT_INPUT_VALUES == 1))
    PRINTF("\nPrinting all Input values other than array contents...\n");
    PRINTF("m = %d\n", m);
    PRINTF("n = %d\n", n);
    PRINTF("lda = %d\n", lda);
    PRINTF("Size of A array (lda*n) = %d\n", lda * n);
    PRINTF("Size of IPIV array (min(m,n)) = %d\n", min(m, n));
  #endif
  
  #if (defined(PRINT_ARRAYS) && (PRINT_ARRAYS == 1))
  // Array to store array name to print.
  char arrayname[20] = "";
  integer arraysize = sizeof(arrayname);
  #endif
  
  #if (defined(PRINT_ARRAYS) && (PRINT_ARRAYS == 1) && \
      defined(PRINT_INPUT_ARRAYS) && (PRINT_INPUT_ARRAYS == 1))
    // Print all input arrays if PRINT_INPUT_ARRAYS macro is enabled
    PRINTF("Printing all Input arrays contents...\n");
    
    // Prints A array contents
    strncpy(arrayname, "A input", arraysize);
    print_array<T>(arrayname, abuff, lda * n);
    strncpy(arrayname, "A ref input", arraysize);
    print_array<T>(arrayname, arefbuff, lda * n);
    
    // Prints IPIV array contents
    strncpy(arrayname, "IPIV input", arraysize);
    print_array<integer>(arrayname, ipivbuff, min(m, n));
    strncpy(arrayname, "IPIV ref input", arraysize);
    print_array<integer>(arrayname, ipivrefbuff, min(m, n));
  #endif
  
  // Call CPP function of GETRF()
  if ((abuff != NULL) && (ipivbuff != NULL)) {
    libflame::getrf<T>(&m, &n, abuff, &lda, ipivbuff, info_cpp);
  }
  
  // Call C function only if input reference buffers are not NULL.
  if ((arefbuff != NULL) && (ipivrefbuff != NULL)) {
    // Call C function of GETRF()
    if (typeid(T) == typeid(float)) {
      getrf_ref = (fptr_NL_LAPACK_getrf)dlsym(lapackModule, "sgetrf_");
    } else if (typeid(T) == typeid(double)) {
      getrf_ref = (fptr_NL_LAPACK_getrf)dlsym(lapackModule, "dgetrf_");
    } else if (typeid(T) == typeid(scomplex)) {
      getrf_ref = (fptr_NL_LAPACK_getrf)dlsym(lapackModule, "cgetrf_");
    } else if (typeid(T) == typeid(dcomplex)) {
      getrf_ref = (fptr_NL_LAPACK_getrf)dlsym(lapackModule, "zgetrf_");
    } else {
      PRINTF("Invalid typename is passed to %s() function template.\n",
             __FUNCTION__);
    }
    
    if (getrf_ref == NULL) {
      PRINTF("Could not get the symbol. Returning...\n");
      return;
    }
    
    getrf_ref(&m, &n, arefbuff, &lda, ipivrefbuff, info_ref);
  }
  
  #if (defined(PRINT_ARRAYS) && (PRINT_ARRAYS == 1) && \
      defined(PRINT_OUTPUT_ARRAYS) && (PRINT_OUTPUT_ARRAYS == 1))
    // Print all output arrays if PRINT_OUTPUT_ARRAYS macro is enabled
    PRINTF("\nPrinting all output arrays contents...\n");
    
    // Prints A array contents
    strncpy(arrayname, "A output", arraysize);
    print_array<T>(arrayname, abuff, lda * n);
    strncpy(arrayname, "A ref output", arraysize);
    print_array<T>(arrayname, arefbuff, lda * n);
    
    // Prints IPIV array contents
    strncpy(arrayname, "IPIV output", arraysize);
    print_array<integer>(arrayname, ipivbuff, min(m, n));
    strncpy(arrayname, "IPIV ref output", arraysize);
    print_array<integer>(arrayname, ipivrefbuff, min(m, n));
  #endif
  PRINTF("%s() Exit...\n", __FUNCTION__);
}

/*! @brief  potrf_internal() function template calls C and CPP based 
      POTRF APIs with valid test values and returns the outout buffers
      and info.
      T can be float, double, scomplex, dcomplex.
 * @details
 * \b Purpose:
    \verbatim
    potrf_internal is function template for potrf() functions.
    T can be float, double, scomplex, dcomplex.
    
    potrf_internal() function template calls C and CPP based POTRF APIs with
    valid test values and returns the outout buffers and info.
    abuff, arefbuff, info_cpp, info_ref will be updated as output.
    If abuff buffer is NULL, then respective CPP function call
    is not performed.
    If arefbuff buffers is NULL, then respective C function call
    is not performed.
    NOTE: It is caller's responsibility to free the buffers after the function
    call.
    \endverbatim
  
 * @params[in] uplo
          UPLO is CHARACTER*1
          = 'U':  Upper triangle of A is stored;
          = 'L':  Lower triangle of A is stored.
 * @params[in] n
          N is INTEGER
          The number of columns of the matrix A.  N >= 0.
 * @params[in] abuff
          A is REAL array, dimension (LDA,N)
          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
          N-by-N upper triangular part of A contains the upper
          triangular part of the matrix A, and the strictly lower
          triangular part of A is not referenced.  If UPLO = 'L', the
          leading N-by-N lower triangular part of A contains the lower
          triangular part of the matrix A, and the strictly upper
          triangular part of A is not referenced.

          On exit, if INFO = 0, the factor U or L from the Cholesky
          factorization A = U**T*U or A = L*L**T.
 * @params[in] lda
          LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,N).
 * @params[in] arefbuff
          arefbuff is REAL array, dimension (LDA,N)
          arefbuff is similar to abuff & used for reference C function.
 * @params[out] info_cpp
          info_cpp is integer to get INFO returned by CPP function.
 * @params[out] info_ref
          info_ref is integer to get INFO returned by reference C function.

 * @return VOID
           Nothing.
 * */
template<typename T>
void potrf_internal(char uplo, integer n, T* abuff, integer lda,
        T* arefbuff, integer* info_cpp, integer* info_ref)
{
  typedef integer (*fptr_NL_LAPACK_potrf)(char* uplo,
                      integer* n, T* a, integer* lda, integer* info);
  fptr_NL_LAPACK_potrf potrf_ref = NULL;
  
  PRINTF("%s() Entry...\n", __FUNCTION__);
  if ((n <= 0) || \
      (lda <= 0) || \
      (info_cpp == NULL) || \
      (info_ref == NULL)) {
    PRINTF("Please pass the valid inputs/buffers. Returning from function.\n");
    return;
  }
  
  if ((abuff == NULL) && (arefbuff == NULL)) {
    PRINTF("Please pass the valid buffers. Returning from function.\n");
    return;
  }
  
  // Print input values other than arrays.
  #if (defined(PRINT_INPUT_VALUES) && (PRINT_INPUT_VALUES == 1))
    PRINTF("\nPrinting all Input values other than array contents...\n");
    PRINTF("uplo = %c\n", uplo);
    PRINTF("n = %d\n", n);
    PRINTF("lda = %d\n", lda);
    PRINTF("Size of A array (lda*n) = %d\n", lda * n);
  #endif
  
  #if (defined(PRINT_ARRAYS) && (PRINT_ARRAYS == 1))
  // Array to store array name to print.
  char arrayname[20] = "";
  integer arraysize = sizeof(arrayname);
  #endif
  
  #if (defined(PRINT_ARRAYS) && (PRINT_ARRAYS == 1) && \
      defined(PRINT_INPUT_ARRAYS) && (PRINT_INPUT_ARRAYS == 1))
    // Print all input arrays if PRINT_INPUT_ARRAYS macro is enabled
    PRINTF("Printing all Input arrays contents...\n");
    
    // Prints A array contents
    strncpy(arrayname, "A input", arraysize);
    print_array<T>(arrayname, abuff, lda * n);
    strncpy(arrayname, "A ref input", arraysize);
    print_array<T>(arrayname, arefbuff, lda * n);
  #endif
  
  // Call CPP function of POTRF()
  if (abuff != NULL) {
    libflame::potrf<T>(&uplo, &n, abuff, &lda, info_cpp);
  }
  
  // Call C function only if input reference buffer is not NULL.
  if (arefbuff != NULL) {
    // Call C function of POTRF()
    if (typeid(T) == typeid(float)) {
      potrf_ref = (fptr_NL_LAPACK_potrf)dlsym(lapackModule, "spotrf_");
    } else if (typeid(T) == typeid(double)) {
      potrf_ref = (fptr_NL_LAPACK_potrf)dlsym(lapackModule, "dpotrf_");
    } else if (typeid(T) == typeid(scomplex)) {
      potrf_ref = (fptr_NL_LAPACK_potrf)dlsym(lapackModule, "cpotrf_");
    } else if (typeid(T) == typeid(dcomplex)) {
      potrf_ref = (fptr_NL_LAPACK_potrf)dlsym(lapackModule, "zpotrf_");
    } else {
      PRINTF("Invalid typename is passed to %s() function template.\n",
             __FUNCTION__);
    }
    
    if (potrf_ref == NULL) {
      PRINTF("Could not get the symbol. Returning...\n");
      return;
    }
    
    potrf_ref(&uplo, &n, arefbuff, &lda, info_ref);
  }
  
  #if (defined(PRINT_ARRAYS) && (PRINT_ARRAYS == 1) && \
      defined(PRINT_OUTPUT_ARRAYS) && (PRINT_OUTPUT_ARRAYS == 1))
    // Print all output arrays if PRINT_OUTPUT_ARRAYS macro is enabled
    PRINTF("\nPrinting all output arrays contents...\n");
    
    // Prints A array contents
    strncpy(arrayname, "A output", arraysize);
    print_array<T>(arrayname, abuff, lda * n);
    strncpy(arrayname, "A ref output", arraysize);
    print_array<T>(arrayname, arefbuff, lda * n);
  #endif
  PRINTF("%s() Exit...\n", __FUNCTION__);
}

/* This function is used to clock the execution time of APIs*/
double fla_test_clock()
{
  double the_time, norm_sec;
  struct timespec ts;

  PRINTF("%s() Entry...\n", __FUNCTION__);
  clock_gettime(CLOCK_MONOTONIC, &ts);

  if (ref_time_sec == 0.0) {
    ref_time_sec = (double) ts.tv_sec;
  }

  norm_sec = (double) ts.tv_sec - ref_time_sec;

  the_time = norm_sec + ts.tv_nsec * 1.0e-9;
  PRINTF("the_time=%lf\n", the_time);
  PRINTF("%s() Exit...\n", __FUNCTION__);
  return the_time;
}
#endif        //  #ifndef LIBFLAME_TEST_HH