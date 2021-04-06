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

void*     bl1_vallocv( uinteger n_elem, uinteger elem_size );
integer*      bl1_iallocv( uinteger n_elem );
float*    bl1_sallocv( uinteger n_elem );
double*   bl1_dallocv( uinteger n_elem );
scomplex* bl1_callocv( uinteger n_elem );
dcomplex* bl1_zallocv( uinteger n_elem );

// --- allocm ---

void*     bl1_vallocm( uinteger m, uinteger n, uinteger elem_size );
integer*      bl1_iallocm( uinteger m, uinteger n );
float*    bl1_sallocm( uinteger m, uinteger n );
double*   bl1_dallocm( uinteger m, uinteger n );
scomplex* bl1_callocm( uinteger m, uinteger n );
dcomplex* bl1_zallocm( uinteger m, uinteger n );

// --- apdiagmv ---

void bl1_sapdiagmv( side1_t side, conj1_t conj, integer m, integer n, float*    x, integer incx, float*    a, integer a_rs, integer a_cs );
void bl1_dapdiagmv( side1_t side, conj1_t conj, integer m, integer n, double*   x, integer incx, double*   a, integer a_rs, integer a_cs );
void bl1_csapdiagmv( side1_t side, conj1_t conj, integer m, integer n, float*    x, integer incx, scomplex* a, integer a_rs, integer a_cs );
void bl1_capdiagmv( side1_t side, conj1_t conj, integer m, integer n, scomplex* x, integer incx, scomplex* a, integer a_rs, integer a_cs );
void bl1_zdapdiagmv( side1_t side, conj1_t conj, integer m, integer n, double*   x, integer incx, dcomplex* a, integer a_rs, integer a_cs );
void bl1_zapdiagmv( side1_t side, conj1_t conj, integer m, integer n, dcomplex* x, integer incx, dcomplex* a, integer a_rs, integer a_cs );

// --- create_contigm ---

void bl1_screate_contigm( integer m, integer n, float*    a_save, integer a_rs_save, integer a_cs_save, float**    a, integer* a_rs, integer* a_cs );
void bl1_dcreate_contigm( integer m, integer n, double*   a_save, integer a_rs_save, integer a_cs_save, double**   a, integer* a_rs, integer* a_cs );
void bl1_ccreate_contigm( integer m, integer n, scomplex* a_save, integer a_rs_save, integer a_cs_save, scomplex** a, integer* a_rs, integer* a_cs );
void bl1_zcreate_contigm( integer m, integer n, dcomplex* a_save, integer a_rs_save, integer a_cs_save, dcomplex** a, integer* a_rs, integer* a_cs );

// --- create_contigmt ---

void bl1_screate_contigmt( trans1_t trans_dims, integer m, integer n, float*    a_save, integer a_rs_save, integer a_cs_save, float**    a, integer* a_rs, integer* a_cs );
void bl1_dcreate_contigmt( trans1_t trans_dims, integer m, integer n, double*   a_save, integer a_rs_save, integer a_cs_save, double**   a, integer* a_rs, integer* a_cs );
void bl1_ccreate_contigmt( trans1_t trans_dims, integer m, integer n, scomplex* a_save, integer a_rs_save, integer a_cs_save, scomplex** a, integer* a_rs, integer* a_cs );
void bl1_zcreate_contigmt( trans1_t trans_dims, integer m, integer n, dcomplex* a_save, integer a_rs_save, integer a_cs_save, dcomplex** a, integer* a_rs, integer* a_cs );

// --- create_contigmr ---

void bl1_screate_contigmr( uplo1_t uplo, integer m, integer n, float*    a_save, integer a_rs_save, integer a_cs_save, float**    a, integer* a_rs, integer* a_cs );
void bl1_dcreate_contigmr( uplo1_t uplo, integer m, integer n, double*   a_save, integer a_rs_save, integer a_cs_save, double**   a, integer* a_rs, integer* a_cs );
void bl1_ccreate_contigmr( uplo1_t uplo, integer m, integer n, scomplex* a_save, integer a_rs_save, integer a_cs_save, scomplex** a, integer* a_rs, integer* a_cs );
void bl1_zcreate_contigmr( uplo1_t uplo, integer m, integer n, dcomplex* a_save, integer a_rs_save, integer a_cs_save, dcomplex** a, integer* a_rs, integer* a_cs );

// --- create_contigmsr ---

void bl1_screate_contigmsr( side1_t side, uplo1_t uplo, integer m, integer n, float*    a_save, integer a_rs_save, integer a_cs_save, float**    a, integer* a_rs, integer* a_cs );
void bl1_dcreate_contigmsr( side1_t side, uplo1_t uplo, integer m, integer n, double*   a_save, integer a_rs_save, integer a_cs_save, double**   a, integer* a_rs, integer* a_cs );
void bl1_ccreate_contigmsr( side1_t side, uplo1_t uplo, integer m, integer n, scomplex* a_save, integer a_rs_save, integer a_cs_save, scomplex** a, integer* a_rs, integer* a_cs );
void bl1_zcreate_contigmsr( side1_t side, uplo1_t uplo, integer m, integer n, dcomplex* a_save, integer a_rs_save, integer a_cs_save, dcomplex** a, integer* a_rs, integer* a_cs );

// --- free_contigm ---

void bl1_sfree_contigm( float*    a_save, integer a_rs_save, integer a_cs_save, float**    a, integer* a_rs, integer* a_cs );
void bl1_dfree_contigm( double*   a_save, integer a_rs_save, integer a_cs_save, double**   a, integer* a_rs, integer* a_cs );
void bl1_cfree_contigm( scomplex* a_save, integer a_rs_save, integer a_cs_save, scomplex** a, integer* a_rs, integer* a_cs );
void bl1_zfree_contigm( dcomplex* a_save, integer a_rs_save, integer a_cs_save, dcomplex** a, integer* a_rs, integer* a_cs );

// --- free_saved_contigm ---

void bl1_sfree_saved_contigm( integer m, integer n, float*    a_save, integer a_rs_save, integer a_cs_save, float**    a, integer* a_rs, integer* a_cs );
void bl1_dfree_saved_contigm( integer m, integer n, double*   a_save, integer a_rs_save, integer a_cs_save, double**   a, integer* a_rs, integer* a_cs );
void bl1_cfree_saved_contigm( integer m, integer n, scomplex* a_save, integer a_rs_save, integer a_cs_save, scomplex** a, integer* a_rs, integer* a_cs );
void bl1_zfree_saved_contigm( integer m, integer n, dcomplex* a_save, integer a_rs_save, integer a_cs_save, dcomplex** a, integer* a_rs, integer* a_cs );

// --- free_saved_contigmr ---

void bl1_sfree_saved_contigmr( uplo1_t uplo, integer m, integer n, float*    a_save, integer a_rs_save, integer a_cs_save, float**    a, integer* a_rs, integer* a_cs );
void bl1_dfree_saved_contigmr( uplo1_t uplo, integer m, integer n, double*   a_save, integer a_rs_save, integer a_cs_save, double**   a, integer* a_rs, integer* a_cs );
void bl1_cfree_saved_contigmr( uplo1_t uplo, integer m, integer n, scomplex* a_save, integer a_rs_save, integer a_cs_save, scomplex** a, integer* a_rs, integer* a_cs );
void bl1_zfree_saved_contigmr( uplo1_t uplo, integer m, integer n, dcomplex* a_save, integer a_rs_save, integer a_cs_save, dcomplex** a, integer* a_rs, integer* a_cs );

// --- free_saved_contigmsr ---

void bl1_sfree_saved_contigmsr( side1_t side, uplo1_t uplo, integer m, integer n, float*    a_save, integer a_rs_save, integer a_cs_save, float**    a, integer* a_rs, integer* a_cs );
void bl1_dfree_saved_contigmsr( side1_t side, uplo1_t uplo, integer m, integer n, double*   a_save, integer a_rs_save, integer a_cs_save, double**   a, integer* a_rs, integer* a_cs );
void bl1_cfree_saved_contigmsr( side1_t side, uplo1_t uplo, integer m, integer n, scomplex* a_save, integer a_rs_save, integer a_cs_save, scomplex** a, integer* a_rs, integer* a_cs );
void bl1_zfree_saved_contigmsr( side1_t side, uplo1_t uplo, integer m, integer n, dcomplex* a_save, integer a_rs_save, integer a_cs_save, dcomplex** a, integer* a_rs, integer* a_cs );

// --- ewinvscalv ---

void bl1_sewinvscalv( conj1_t conj, integer n, float*    x, integer incx, float*    y, integer incy );
void bl1_dewinvscalv( conj1_t conj, integer n, double*   x, integer incx, double*   y, integer incy );
void bl1_csewinvscalv( conj1_t conj, integer n, float*    x, integer incx, scomplex* y, integer incy );
void bl1_cewinvscalv( conj1_t conj, integer n, scomplex* x, integer incx, scomplex* y, integer incy );
void bl1_zdewinvscalv( conj1_t conj, integer n, double*   x, integer incx, dcomplex* y, integer incy );
void bl1_zewinvscalv( conj1_t conj, integer n, dcomplex* x, integer incx, dcomplex* y, integer incy );

// --- ewscalmt ---

void bl1_sewinvscalmt( trans1_t trans, integer m, integer n, float*    a, integer a_rs, integer a_cs, float*    b, integer b_rs, integer b_cs );
void bl1_dewinvscalmt( trans1_t trans, integer m, integer n, double*   a, integer a_rs, integer a_cs, double*   b, integer b_rs, integer b_cs );
void bl1_csewinvscalmt( trans1_t trans, integer m, integer n, float*    a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs );
void bl1_cewinvscalmt( trans1_t trans, integer m, integer n, scomplex* a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs );
void bl1_zdewinvscalmt( trans1_t trans, integer m, integer n, double*   a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs );
void bl1_zewinvscalmt( trans1_t trans, integer m, integer n, dcomplex* a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs );

// --- ewscalv ---

void bl1_sewscalv( conj1_t conj, integer n, float*    x, integer incx, float*    y, integer incy );
void bl1_dewscalv( conj1_t conj, integer n, double*   x, integer incx, double*   y, integer incy );
void bl1_csewscalv( conj1_t conj, integer n, float*    x, integer incx, scomplex* y, integer incy );
void bl1_cewscalv( conj1_t conj, integer n, scomplex* x, integer incx, scomplex* y, integer incy );
void bl1_zdewscalv( conj1_t conj, integer n, double*   x, integer incx, dcomplex* y, integer incy );
void bl1_zewscalv( conj1_t conj, integer n, dcomplex* x, integer incx, dcomplex* y, integer incy );

// --- ewscalmt ---

void bl1_sewscalmt( trans1_t trans, integer m, integer n, float*    a, integer a_rs, integer a_cs, float*    b, integer b_rs, integer b_cs );
void bl1_dewscalmt( trans1_t trans, integer m, integer n, double*   a, integer a_rs, integer a_cs, double*   b, integer b_rs, integer b_cs );
void bl1_csewscalmt( trans1_t trans, integer m, integer n, float*    a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs );
void bl1_cewscalmt( trans1_t trans, integer m, integer n, scomplex* a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs );
void bl1_zdewscalmt( trans1_t trans, integer m, integer n, double*   a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs );
void bl1_zewscalmt( trans1_t trans, integer m, integer n, dcomplex* a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs );

// --- free ---

void bl1_vfree( void*     p );
void bl1_ifree( integer*      p );
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

void bl1_sinvertv( conj1_t conj, integer n, float*    x, integer incx );
void bl1_dinvertv( conj1_t conj, integer n, double*   x, integer incx );
void bl1_cinvertv( conj1_t conj, integer n, scomplex* x, integer incx );
void bl1_zinvertv( conj1_t conj, integer n, dcomplex* x, integer incx );

// --- ident ---

void bl1_sident( integer m, float*    a, integer a_rs, integer a_cs );
void bl1_dident( integer m, double*   a, integer a_rs, integer a_cs );
void bl1_cident( integer m, scomplex* a, integer a_rs, integer a_cs );
void bl1_zident( integer m, dcomplex* a, integer a_rs, integer a_cs );

// --- maxabsv ---

void bl1_smaxabsv( integer n, float*    x, integer incx, float*  maxabs );
void bl1_dmaxabsv( integer n, double*   x, integer incx, double* maxabs );
void bl1_cmaxabsv( integer n, scomplex* x, integer incx, float*  maxabs );
void bl1_zmaxabsv( integer n, dcomplex* x, integer incx, double* maxabs );

// --- maxabsm ---

void bl1_smaxabsm( integer m, integer n, float*    a, integer a_rs, integer a_cs, float*  maxabs );
void bl1_dmaxabsm( integer m, integer n, double*   a, integer a_rs, integer a_cs, double* maxabs );
void bl1_cmaxabsm( integer m, integer n, scomplex* a, integer a_rs, integer a_cs, float*  maxabs );
void bl1_zmaxabsm( integer m, integer n, dcomplex* a, integer a_rs, integer a_cs, double* maxabs );

// --- maxabsmr ---

void bl1_smaxabsmr( uplo1_t uplo, integer m, integer n, float*    a, integer a_rs, integer a_cs, float*  maxabs );
void bl1_dmaxabsmr( uplo1_t uplo, integer m, integer n, double*   a, integer a_rs, integer a_cs, double* maxabs );
void bl1_cmaxabsmr( uplo1_t uplo, integer m, integer n, scomplex* a, integer a_rs, integer a_cs, float*  maxabs );
void bl1_zmaxabsmr( uplo1_t uplo, integer m, integer n, dcomplex* a, integer a_rs, integer a_cs, double* maxabs );

// --- rands ---

void bl1_srands( float*    alpha );
void bl1_drands( double*   alpha );
void bl1_crands( scomplex* alpha );
void bl1_zrands( dcomplex* alpha );

// --- randv ---

void bl1_srandv( integer n, float*    x, integer incx );
void bl1_drandv( integer n, double*   x, integer incx );
void bl1_crandv( integer n, scomplex* x, integer incx );
void bl1_zrandv( integer n, dcomplex* x, integer incx );

// --- randm ---

void bl1_srandm( integer m, integer n, float*    a, integer a_rs, integer a_cs );
void bl1_drandm( integer m, integer n, double*   a, integer a_rs, integer a_cs );
void bl1_crandm( integer m, integer n, scomplex* a, integer a_rs, integer a_cs );
void bl1_zrandm( integer m, integer n, dcomplex* a, integer a_rs, integer a_cs );

// --- randmr ---
void bl1_srandmr( uplo1_t uplo, diag1_t diag, integer m, integer n, float*    a, integer a_rs, integer a_cs );
void bl1_drandmr( uplo1_t uplo, diag1_t diag, integer m, integer n, double*   a, integer a_rs, integer a_cs );
void bl1_crandmr( uplo1_t uplo, diag1_t diag, integer m, integer n, scomplex* a, integer a_rs, integer a_cs );
void bl1_zrandmr( uplo1_t uplo, diag1_t diag, integer m, integer n, dcomplex* a, integer a_rs, integer a_cs );

// --- set_contig_strides ---

void bl1_set_contig_strides( integer m, integer n, integer* rs, integer* cs );

// --- set_dims_with_side ---

void bl1_set_dim_with_side( side1_t side, integer m, integer n, integer* dim_new );

// --- set_dims_with_trans ---

void bl1_set_dims_with_trans( trans1_t trans, integer m, integer n, integer* m_new, integer* n_new );

// --- setv ---

void bl1_isetv( integer m, integer*      sigma, integer*      x, integer incx );
void bl1_ssetv( integer m, float*    sigma, float*    x, integer incx );
void bl1_dsetv( integer m, double*   sigma, double*   x, integer incx );
void bl1_csetv( integer m, scomplex* sigma, scomplex* x, integer incx );
void bl1_zsetv( integer m, dcomplex* sigma, dcomplex* x, integer incx );

// --- setm ---

void bl1_isetm( integer m, integer n, integer*      sigma, integer*      a, integer a_rs, integer a_cs );
void bl1_ssetm( integer m, integer n, float*    sigma, float*    a, integer a_rs, integer a_cs );
void bl1_dsetm( integer m, integer n, double*   sigma, double*   a, integer a_rs, integer a_cs );
void bl1_csetm( integer m, integer n, scomplex* sigma, scomplex* a, integer a_rs, integer a_cs );
void bl1_zsetm( integer m, integer n, dcomplex* sigma, dcomplex* a, integer a_rs, integer a_cs );

// --- setmr ---

void bl1_ssetmr( uplo1_t uplo, integer m, integer n, float*    sigma, float*    a, integer a_rs, integer a_cs );
void bl1_dsetmr( uplo1_t uplo, integer m, integer n, double*   sigma, double*   a, integer a_rs, integer a_cs );
void bl1_csetmr( uplo1_t uplo, integer m, integer n, scomplex* sigma, scomplex* a, integer a_rs, integer a_cs );
void bl1_zsetmr( uplo1_t uplo, integer m, integer n, dcomplex* sigma, dcomplex* a, integer a_rs, integer a_cs );

// --- setdiag ---

void bl1_isetdiag( integer offset, integer m, integer n, integer*      sigma, integer*      a, integer a_rs, integer a_cs );
void bl1_ssetdiag( integer offset, integer m, integer n, float*    sigma, float*    a, integer a_rs, integer a_cs );
void bl1_dsetdiag( integer offset, integer m, integer n, double*   sigma, double*   a, integer a_rs, integer a_cs );
void bl1_csetdiag( integer offset, integer m, integer n, scomplex* sigma, scomplex* a, integer a_rs, integer a_cs );
void bl1_zsetdiag( integer offset, integer m, integer n, dcomplex* sigma, dcomplex* a, integer a_rs, integer a_cs );

// --- scalediag ---

void bl1_sscalediag( conj1_t conj, integer offset, integer m, integer n, float*    sigma, float*    a, integer a_rs, integer a_cs );
void bl1_dscalediag( conj1_t conj, integer offset, integer m, integer n, double*   sigma, double*   a, integer a_rs, integer a_cs );
void bl1_cscalediag( conj1_t conj, integer offset, integer m, integer n, scomplex* sigma, scomplex* a, integer a_rs, integer a_cs );
void bl1_zscalediag( conj1_t conj, integer offset, integer m, integer n, dcomplex* sigma, dcomplex* a, integer a_rs, integer a_cs );
void bl1_csscalediag( conj1_t conj, integer offset, integer m, integer n, float*    sigma, scomplex* a, integer a_rs, integer a_cs );
void bl1_zdscalediag( conj1_t conj, integer offset, integer m, integer n, double*   sigma, dcomplex* a, integer a_rs, integer a_cs );

// --- shiftdiag ---

void bl1_sshiftdiag( conj1_t conj, integer offset, integer m, integer n, float*    sigma, float*    a, integer a_rs, integer a_cs );
void bl1_dshiftdiag( conj1_t conj, integer offset, integer m, integer n, double*   sigma, double*   a, integer a_rs, integer a_cs );
void bl1_cshiftdiag( conj1_t conj, integer offset, integer m, integer n, scomplex* sigma, scomplex* a, integer a_rs, integer a_cs );
void bl1_zshiftdiag( conj1_t conj, integer offset, integer m, integer n, dcomplex* sigma, dcomplex* a, integer a_rs, integer a_cs );
void bl1_csshiftdiag( conj1_t conj, integer offset, integer m, integer n, float*    sigma, scomplex* a, integer a_rs, integer a_cs );
void bl1_zdshiftdiag( conj1_t conj, integer offset, integer m, integer n, double*   sigma, dcomplex* a, integer a_rs, integer a_cs );

// --- symmize ---

void bl1_ssymmize( conj1_t conj, uplo1_t uplo, integer m, float*    a, integer a_rs, integer a_cs );
void bl1_dsymmize( conj1_t conj, uplo1_t uplo, integer m, double*   a, integer a_rs, integer a_cs );
void bl1_csymmize( conj1_t conj, uplo1_t uplo, integer m, scomplex* a, integer a_rs, integer a_cs );
void bl1_zsymmize( conj1_t conj, uplo1_t uplo, integer m, dcomplex* a, integer a_rs, integer a_cs );

