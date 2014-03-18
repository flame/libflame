
// --- Abort prototypes --------------------------------------------------------

void bl1_abort( void );
void bl1_abort_msg( char* message );

// --- Parameter-mapping prototypes --------------------------------------------

void bl1_param_map_to_netlib_trans( trans1_t blis_trans, void* blas_trans );
void bl1_param_map_to_netlib_uplo(  uplo1_t  blis_uplo,  void* blas_uplo );
void bl1_param_map_to_netlib_side(  side1_t  blis_side,  void* blas_side );
void bl1_param_map_to_netlib_diag(  diag1_t  blis_diag,  void* blas_diag );

