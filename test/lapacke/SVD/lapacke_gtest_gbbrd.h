#ifndef LAPACKE_GTEST_GBBRD_H
#define LAPACKE_GTEST_GBBRD_H

#define gbbrd_free() \
       free (ab   ); \
       free (abref); \
       free (c   ); \
       free (cref); \
       free (e   ); \
       free (eref); \
       free (d   ); \
       free (dref); \
       free (q   ); \
       free (qref); \
       free (pt   ); \
       free (ptref); \
       free (work  ); \
       free (workref)

class gbbrd_double_parameters{
	
   public:
   /* input params to the API **/
   int 	matrix_layout;
   char vect;
   int 	m;
   int 	n;
   int 	ncc;
   int 	kl;
   int 	ku;
   double *ab, *abref;
   int 	ldab;
   double *d, *dref;
   double *e, *eref;
   double *q, *qref;
   int 	ldq;
   double *pt, *ptref;
   int 	ldpt;
   double *c, *cref;
   int 	ldc;
   double *work, *workref;
   int info, inforef;

   public: 
      gbbrd_double_parameters (int matrix_layout_i, char vect_i, int m_i,
                               int n_i, int ncc_i, int kl_i, int ku_i,
                               int ldab_i, int ldq_i, int ldpt_i, int ldc_i );
      ~gbbrd_double_parameters ();

}; /* end of gbbrd_double_parameters  class definition */


/* 'gbbrd_float_parameters' class definition  */

class gbbrd_float_parameters{
	
   public:
   /* input params to the API **/
   int 	matrix_layout;
   char vect;
   int 	m;
   int 	n;
   int 	ncc;
   int 	kl;
   int 	ku;
   float *ab, *abref;
   int 	ldab;
   float *d, *dref;
   float *e, *eref;
   float *q, *qref;
   int 	ldq;
   float *pt, *ptref;
   int 	ldpt;
   float *c, *cref;
   int 	ldc;
   float *work, *workref;
   int info, inforef;

   public: 
      gbbrd_float_parameters (int matrix_layout_i, char vect_i, int m_i,
                               int n_i, int ncc_i, int kl_i, int ku_i,
                               int ldab_i, int ldq_i, int ldpt_i, int ldc_i );
      ~gbbrd_float_parameters ();

}; /* end of gbbrd_float_parameters  class definition */

class gbbrd_scomplex_parameters{
	
   public:
   /* input params to the API **/
   int 	matrix_layout;
   char vect;
   int 	m;
   int 	n;
   int 	ncc;
   int 	kl;
   int 	ku;
   lapack_complex_float *ab, *abref;
   int 	ldab;
   float *d, *dref;
   float *e, *eref;
   lapack_complex_float *q, *qref;
   int 	ldq;
   lapack_complex_float *pt, *ptref;
   int 	ldpt;
   lapack_complex_float *c, *cref;
   int 	ldc;
   lapack_complex_float *work, *workref;
   int info, inforef;

   public: 
      gbbrd_scomplex_parameters (int matrix_layout_i, char vect_i, int m_i,
                               int n_i, int ncc_i, int kl_i, int ku_i,
                               int ldab_i, int ldq_i, int ldpt_i, int ldc_i );
      ~gbbrd_scomplex_parameters ();

}; /* end of gbbrd_scomplex_parameters  class definition */

class gbbrd_dcomplex_parameters{
	
   public:
   /* input params to the API **/
   int 	matrix_layout;
   char vect;
   int 	m;
   int 	n;
   int 	ncc;
   int 	kl;
   int 	ku;
   lapack_complex_double *ab, *abref;
   int 	ldab;
   double *d, *dref;
   double *e, *eref;
   lapack_complex_double *q, *qref;
   int 	ldq;
   lapack_complex_double *pt, *ptref;
   int 	ldpt;
   lapack_complex_double *c, *cref;
   int 	ldc;
   lapack_complex_double *work, *workref;
   int info, inforef;

   public: 
      gbbrd_dcomplex_parameters (int matrix_layout_i, char vect_i, int m_i,
                               int n_i, int ncc_i, int kl_i, int ku_i,
                               int ldab_i, int ldq_i, int ldpt_i, int ldc_i );
      ~gbbrd_dcomplex_parameters ();

}; /* end of gbbrd_dcomplex_parameters  class definition */


#endif