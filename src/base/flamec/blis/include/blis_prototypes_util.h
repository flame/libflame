/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// --- Utility-level BLAS-like prototypes --------------------------------------

// --- constant-generating functions ---

float    bl1_s2( void );
double   bl1_d2( void );
scomplex bl1_c2( void );
dcomplex bl1_z2( void );
float    bl1_s1( void );
double   bl1_d1( void );
scomplex bl1_c1( void );
dcomplex bl1_z1( void );
float    bl1_s1h( void );
double   bl1_d1h( void );
scomplex bl1_c1h( void );
dcomplex bl1_z1h( void );
float    bl1_s0( void );
double   bl1_d0( void );
scomplex bl1_c0( void );
dcomplex bl1_z0( void );
float    bl1_sm1h( void );
double   bl1_dm1h( void );
scomplex bl1_cm1h( void );
dcomplex bl1_zm1h( void );
float    bl1_sm1( void );
double   bl1_dm1( void );
scomplex bl1_cm1( void );
dcomplex bl1_zm1( void );
float    bl1_sm2( void );
double   bl1_dm2( void );
scomplex bl1_cm2( void );
dcomplex bl1_zm2( void );

// --- allocv ---

void*     bl1_vallocv( unsigned int n_elem, unsigned int elem_size );
int*      bl1_iallocv( unsigned int n_elem );
float*    bl1_sallocv( unsigned int n_elem );
double*   bl1_dallocv( unsigned int n_elem );
scomplex* bl1_callocv( unsigned int n_elem );
dcomplex* bl1_zallocv( unsigned int n_elem );

// --- allocm ---

void*     bl1_vallocm( unsigned int m, unsigned int n, unsigned int elem_size );
int*      bl1_iallocm( unsigned int m, unsigned int n );
float*    bl1_sallocm( unsigned int m, unsigned int n );
double*   bl1_dallocm( unsigned int m, unsigned int n );
scomplex* bl1_callocm( unsigned int m, unsigned int n );
dcomplex* bl1_zallocm( unsigned int m, unsigned int n );

// --- apdiagmv ---

void bl1_sapdiagmv( side1_t side, conj1_t conj, int m, int n, float*    x, int incx, float*    a, int a_rs, int a_cs );
void bl1_dapdiagmv( side1_t side, conj1_t conj, int m, int n, double*   x, int incx, double*   a, int a_rs, int a_cs );
void bl1_csapdiagmv( side1_t side, conj1_t conj, int m, int n, float*    x, int incx, scomplex* a, int a_rs, int a_cs );
void bl1_capdiagmv( side1_t side, conj1_t conj, int m, int n, scomplex* x, int incx, scomplex* a, int a_rs, int a_cs );
void bl1_zdapdiagmv( side1_t side, conj1_t conj, int m, int n, double*   x, int incx, dcomplex* a, int a_rs, int a_cs );
void bl1_zapdiagmv( side1_t side, conj1_t conj, int m, int n, dcomplex* x, int incx, dcomplex* a, int a_rs, int a_cs );

// --- create_contigm ---

void bl1_screate_contigm( int m, int n, float*    a_save, int a_rs_save, int a_cs_save, float**    a, int* a_rs, int* a_cs );
void bl1_dcreate_contigm( int m, int n, double*   a_save, int a_rs_save, int a_cs_save, double**   a, int* a_rs, int* a_cs );
void bl1_ccreate_contigm( int m, int n, scomplex* a_save, int a_rs_save, int a_cs_save, scomplex** a, int* a_rs, int* a_cs );
void bl1_zcreate_contigm( int m, int n, dcomplex* a_save, int a_rs_save, int a_cs_save, dcomplex** a, int* a_rs, int* a_cs );

// --- create_contigmt ---

void bl1_screate_contigmt( trans1_t trans_dims, int m, int n, float*    a_save, int a_rs_save, int a_cs_save, float**    a, int* a_rs, int* a_cs );
void bl1_dcreate_contigmt( trans1_t trans_dims, int m, int n, double*   a_save, int a_rs_save, int a_cs_save, double**   a, int* a_rs, int* a_cs );
void bl1_ccreate_contigmt( trans1_t trans_dims, int m, int n, scomplex* a_save, int a_rs_save, int a_cs_save, scomplex** a, int* a_rs, int* a_cs );
void bl1_zcreate_contigmt( trans1_t trans_dims, int m, int n, dcomplex* a_save, int a_rs_save, int a_cs_save, dcomplex** a, int* a_rs, int* a_cs );

// --- create_contigmr ---

void bl1_screate_contigmr( uplo1_t uplo, int m, int n, float*    a_save, int a_rs_save, int a_cs_save, float**    a, int* a_rs, int* a_cs );
void bl1_dcreate_contigmr( uplo1_t uplo, int m, int n, double*   a_save, int a_rs_save, int a_cs_save, double**   a, int* a_rs, int* a_cs );
void bl1_ccreate_contigmr( uplo1_t uplo, int m, int n, scomplex* a_save, int a_rs_save, int a_cs_save, scomplex** a, int* a_rs, int* a_cs );
void bl1_zcreate_contigmr( uplo1_t uplo, int m, int n, dcomplex* a_save, int a_rs_save, int a_cs_save, dcomplex** a, int* a_rs, int* a_cs );

// --- create_contigmsr ---

void bl1_screate_contigmsr( side1_t side, uplo1_t uplo, int m, int n, float*    a_save, int a_rs_save, int a_cs_save, float**    a, int* a_rs, int* a_cs );
void bl1_dcreate_contigmsr( side1_t side, uplo1_t uplo, int m, int n, double*   a_save, int a_rs_save, int a_cs_save, double**   a, int* a_rs, int* a_cs );
void bl1_ccreate_contigmsr( side1_t side, uplo1_t uplo, int m, int n, scomplex* a_save, int a_rs_save, int a_cs_save, scomplex** a, int* a_rs, int* a_cs );
void bl1_zcreate_contigmsr( side1_t side, uplo1_t uplo, int m, int n, dcomplex* a_save, int a_rs_save, int a_cs_save, dcomplex** a, int* a_rs, int* a_cs );

// --- free_contigm ---

void bl1_sfree_contigm( float*    a_save, int a_rs_save, int a_cs_save, float**    a, int* a_rs, int* a_cs );
void bl1_dfree_contigm( double*   a_save, int a_rs_save, int a_cs_save, double**   a, int* a_rs, int* a_cs );
void bl1_cfree_contigm( scomplex* a_save, int a_rs_save, int a_cs_save, scomplex** a, int* a_rs, int* a_cs );
void bl1_zfree_contigm( dcomplex* a_save, int a_rs_save, int a_cs_save, dcomplex** a, int* a_rs, int* a_cs );

// --- free_saved_contigm ---

void bl1_sfree_saved_contigm( int m, int n, float*    a_save, int a_rs_save, int a_cs_save, float**    a, int* a_rs, int* a_cs );
void bl1_dfree_saved_contigm( int m, int n, double*   a_save, int a_rs_save, int a_cs_save, double**   a, int* a_rs, int* a_cs );
void bl1_cfree_saved_contigm( int m, int n, scomplex* a_save, int a_rs_save, int a_cs_save, scomplex** a, int* a_rs, int* a_cs );
void bl1_zfree_saved_contigm( int m, int n, dcomplex* a_save, int a_rs_save, int a_cs_save, dcomplex** a, int* a_rs, int* a_cs );

// --- free_saved_contigmr ---

void bl1_sfree_saved_contigmr( uplo1_t uplo, int m, int n, float*    a_save, int a_rs_save, int a_cs_save, float**    a, int* a_rs, int* a_cs );
void bl1_dfree_saved_contigmr( uplo1_t uplo, int m, int n, double*   a_save, int a_rs_save, int a_cs_save, double**   a, int* a_rs, int* a_cs );
void bl1_cfree_saved_contigmr( uplo1_t uplo, int m, int n, scomplex* a_save, int a_rs_save, int a_cs_save, scomplex** a, int* a_rs, int* a_cs );
void bl1_zfree_saved_contigmr( uplo1_t uplo, int m, int n, dcomplex* a_save, int a_rs_save, int a_cs_save, dcomplex** a, int* a_rs, int* a_cs );

// --- free_saved_contigmsr ---

void bl1_sfree_saved_contigmsr( side1_t side, uplo1_t uplo, int m, int n, float*    a_save, int a_rs_save, int a_cs_save, float**    a, int* a_rs, int* a_cs );
void bl1_dfree_saved_contigmsr( side1_t side, uplo1_t uplo, int m, int n, double*   a_save, int a_rs_save, int a_cs_save, double**   a, int* a_rs, int* a_cs );
void bl1_cfree_saved_contigmsr( side1_t side, uplo1_t uplo, int m, int n, scomplex* a_save, int a_rs_save, int a_cs_save, scomplex** a, int* a_rs, int* a_cs );
void bl1_zfree_saved_contigmsr( side1_t side, uplo1_t uplo, int m, int n, dcomplex* a_save, int a_rs_save, int a_cs_save, dcomplex** a, int* a_rs, int* a_cs );

// --- ewinvscalv ---

void bl1_sewinvscalv( conj1_t conj, int n, float*    x, int incx, float*    y, int incy );
void bl1_dewinvscalv( conj1_t conj, int n, double*   x, int incx, double*   y, int incy );
void bl1_csewinvscalv( conj1_t conj, int n, float*    x, int incx, scomplex* y, int incy );
void bl1_cewinvscalv( conj1_t conj, int n, scomplex* x, int incx, scomplex* y, int incy );
void bl1_zdewinvscalv( conj1_t conj, int n, double*   x, int incx, dcomplex* y, int incy );
void bl1_zewinvscalv( conj1_t conj, int n, dcomplex* x, int incx, dcomplex* y, int incy );

// --- ewscalmt ---

void bl1_sewinvscalmt( trans1_t trans, int m, int n, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_dewinvscalmt( trans1_t trans, int m, int n, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_csewinvscalmt( trans1_t trans, int m, int n, float*    a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_cewinvscalmt( trans1_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_zdewinvscalmt( trans1_t trans, int m, int n, double*   a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );
void bl1_zewinvscalmt( trans1_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

// --- ewscalv ---

void bl1_sewscalv( conj1_t conj, int n, float*    x, int incx, float*    y, int incy );
void bl1_dewscalv( conj1_t conj, int n, double*   x, int incx, double*   y, int incy );
void bl1_csewscalv( conj1_t conj, int n, float*    x, int incx, scomplex* y, int incy );
void bl1_cewscalv( conj1_t conj, int n, scomplex* x, int incx, scomplex* y, int incy );
void bl1_zdewscalv( conj1_t conj, int n, double*   x, int incx, dcomplex* y, int incy );
void bl1_zewscalv( conj1_t conj, int n, dcomplex* x, int incx, dcomplex* y, int incy );

// --- ewscalmt ---

void bl1_sewscalmt( trans1_t trans, int m, int n, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_dewscalmt( trans1_t trans, int m, int n, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_csewscalmt( trans1_t trans, int m, int n, float*    a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_cewscalmt( trans1_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_zdewscalmt( trans1_t trans, int m, int n, double*   a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );
void bl1_zewscalmt( trans1_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

// --- free ---

void bl1_vfree( void*     p );
void bl1_ifree( int*      p );
void bl1_sfree( float*    p );
void bl1_dfree( double*   p );
void bl1_cfree( scomplex* p );
void bl1_zfree( dcomplex* p );

// --- inverts ---

void bl1_sinverts( conj1_t conj, float*    alpha );
void bl1_dinverts( conj1_t conj, double*   alpha );
void bl1_cinverts( conj1_t conj, scomplex* alpha );
void bl1_zinverts( conj1_t conj, dcomplex* alpha );

// --- invert2s ---

void bl1_sinvert2s( conj1_t conj, float*    alpha, float*    beta );
void bl1_dinvert2s( conj1_t conj, double*   alpha, double*   beta );
void bl1_cinvert2s( conj1_t conj, scomplex* alpha, scomplex* beta );
void bl1_zinvert2s( conj1_t conj, dcomplex* alpha, dcomplex* beta );

// --- invertv ---

void bl1_sinvertv( conj1_t conj, int n, float*    x, int incx );
void bl1_dinvertv( conj1_t conj, int n, double*   x, int incx );
void bl1_cinvertv( conj1_t conj, int n, scomplex* x, int incx );
void bl1_zinvertv( conj1_t conj, int n, dcomplex* x, int incx );

// --- ident ---

void bl1_sident( int m, float*    a, int a_rs, int a_cs );
void bl1_dident( int m, double*   a, int a_rs, int a_cs );
void bl1_cident( int m, scomplex* a, int a_rs, int a_cs );
void bl1_zident( int m, dcomplex* a, int a_rs, int a_cs );

// --- maxabsv ---

void bl1_smaxabsv( int n, float*    x, int incx, float*  maxabs );
void bl1_dmaxabsv( int n, double*   x, int incx, double* maxabs );
void bl1_cmaxabsv( int n, scomplex* x, int incx, float*  maxabs );
void bl1_zmaxabsv( int n, dcomplex* x, int incx, double* maxabs );

// --- maxabsm ---

void bl1_smaxabsm( int m, int n, float*    a, int a_rs, int a_cs, float*  maxabs );
void bl1_dmaxabsm( int m, int n, double*   a, int a_rs, int a_cs, double* maxabs );
void bl1_cmaxabsm( int m, int n, scomplex* a, int a_rs, int a_cs, float*  maxabs );
void bl1_zmaxabsm( int m, int n, dcomplex* a, int a_rs, int a_cs, double* maxabs );

// --- maxabsmr ---

void bl1_smaxabsmr( uplo1_t uplo, int m, int n, float*    a, int a_rs, int a_cs, float*  maxabs );
void bl1_dmaxabsmr( uplo1_t uplo, int m, int n, double*   a, int a_rs, int a_cs, double* maxabs );
void bl1_cmaxabsmr( uplo1_t uplo, int m, int n, scomplex* a, int a_rs, int a_cs, float*  maxabs );
void bl1_zmaxabsmr( uplo1_t uplo, int m, int n, dcomplex* a, int a_rs, int a_cs, double* maxabs );

// --- rands ---

void bl1_srands( float*    alpha );
void bl1_drands( double*   alpha );
void bl1_crands( scomplex* alpha );
void bl1_zrands( dcomplex* alpha );

// --- randv ---

void bl1_srandv( int n, float*    x, int incx );
void bl1_drandv( int n, double*   x, int incx );
void bl1_crandv( int n, scomplex* x, int incx );
void bl1_zrandv( int n, dcomplex* x, int incx );

// --- randm ---

void bl1_srandm( int m, int n, float*    a, int a_rs, int a_cs );
void bl1_drandm( int m, int n, double*   a, int a_rs, int a_cs );
void bl1_crandm( int m, int n, scomplex* a, int a_rs, int a_cs );
void bl1_zrandm( int m, int n, dcomplex* a, int a_rs, int a_cs );

// --- randmr ---
void bl1_srandmr( uplo1_t uplo, diag1_t diag, int m, int n, float*    a, int a_rs, int a_cs );
void bl1_drandmr( uplo1_t uplo, diag1_t diag, int m, int n, double*   a, int a_rs, int a_cs );
void bl1_crandmr( uplo1_t uplo, diag1_t diag, int m, int n, scomplex* a, int a_rs, int a_cs );
void bl1_zrandmr( uplo1_t uplo, diag1_t diag, int m, int n, dcomplex* a, int a_rs, int a_cs );

// --- set_contig_strides ---

void bl1_set_contig_strides( int m, int n, int* rs, int* cs );

// --- set_dims_with_side ---

void bl1_set_dim_with_side( side1_t side, int m, int n, int* dim_new );

// --- set_dims_with_trans ---

void bl1_set_dims_with_trans( trans1_t trans, int m, int n, int* m_new, int* n_new );

// --- setv ---

void bl1_isetv( int m, int*      sigma, int*      x, int incx );
void bl1_ssetv( int m, float*    sigma, float*    x, int incx );
void bl1_dsetv( int m, double*   sigma, double*   x, int incx );
void bl1_csetv( int m, scomplex* sigma, scomplex* x, int incx );
void bl1_zsetv( int m, dcomplex* sigma, dcomplex* x, int incx );

// --- setm ---

void bl1_isetm( int m, int n, int*      sigma, int*      a, int a_rs, int a_cs );
void bl1_ssetm( int m, int n, float*    sigma, float*    a, int a_rs, int a_cs );
void bl1_dsetm( int m, int n, double*   sigma, double*   a, int a_rs, int a_cs );
void bl1_csetm( int m, int n, scomplex* sigma, scomplex* a, int a_rs, int a_cs );
void bl1_zsetm( int m, int n, dcomplex* sigma, dcomplex* a, int a_rs, int a_cs );

// --- setmr ---

void bl1_ssetmr( uplo1_t uplo, int m, int n, float*    sigma, float*    a, int a_rs, int a_cs );
void bl1_dsetmr( uplo1_t uplo, int m, int n, double*   sigma, double*   a, int a_rs, int a_cs );
void bl1_csetmr( uplo1_t uplo, int m, int n, scomplex* sigma, scomplex* a, int a_rs, int a_cs );
void bl1_zsetmr( uplo1_t uplo, int m, int n, dcomplex* sigma, dcomplex* a, int a_rs, int a_cs );

// --- setdiag ---

void bl1_isetdiag( int offset, int m, int n, int*      sigma, int*      a, int a_rs, int a_cs );
void bl1_ssetdiag( int offset, int m, int n, float*    sigma, float*    a, int a_rs, int a_cs );
void bl1_dsetdiag( int offset, int m, int n, double*   sigma, double*   a, int a_rs, int a_cs );
void bl1_csetdiag( int offset, int m, int n, scomplex* sigma, scomplex* a, int a_rs, int a_cs );
void bl1_zsetdiag( int offset, int m, int n, dcomplex* sigma, dcomplex* a, int a_rs, int a_cs );

// --- scalediag ---

void bl1_sscalediag( conj1_t conj, int offset, int m, int n, float*    sigma, float*    a, int a_rs, int a_cs );
void bl1_dscalediag( conj1_t conj, int offset, int m, int n, double*   sigma, double*   a, int a_rs, int a_cs );
void bl1_cscalediag( conj1_t conj, int offset, int m, int n, scomplex* sigma, scomplex* a, int a_rs, int a_cs );
void bl1_zscalediag( conj1_t conj, int offset, int m, int n, dcomplex* sigma, dcomplex* a, int a_rs, int a_cs );
void bl1_csscalediag( conj1_t conj, int offset, int m, int n, float*    sigma, scomplex* a, int a_rs, int a_cs );
void bl1_zdscalediag( conj1_t conj, int offset, int m, int n, double*   sigma, dcomplex* a, int a_rs, int a_cs );

// --- shiftdiag ---

void bl1_sshiftdiag( conj1_t conj, int offset, int m, int n, float*    sigma, float*    a, int a_rs, int a_cs );
void bl1_dshiftdiag( conj1_t conj, int offset, int m, int n, double*   sigma, double*   a, int a_rs, int a_cs );
void bl1_cshiftdiag( conj1_t conj, int offset, int m, int n, scomplex* sigma, scomplex* a, int a_rs, int a_cs );
void bl1_zshiftdiag( conj1_t conj, int offset, int m, int n, dcomplex* sigma, dcomplex* a, int a_rs, int a_cs );
void bl1_csshiftdiag( conj1_t conj, int offset, int m, int n, float*    sigma, scomplex* a, int a_rs, int a_cs );
void bl1_zdshiftdiag( conj1_t conj, int offset, int m, int n, double*   sigma, dcomplex* a, int a_rs, int a_cs );

// --- symmize ---

void bl1_ssymmize( conj1_t conj, uplo1_t uplo, int m, float*    a, int a_rs, int a_cs );
void bl1_dsymmize( conj1_t conj, uplo1_t uplo, int m, double*   a, int a_rs, int a_cs );
void bl1_csymmize( conj1_t conj, uplo1_t uplo, int m, scomplex* a, int a_rs, int a_cs );
void bl1_zsymmize( conj1_t conj, uplo1_t uplo, int m, dcomplex* a, int a_rs, int a_cs );

