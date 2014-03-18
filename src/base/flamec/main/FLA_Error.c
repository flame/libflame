
#include "FLAME.h"

// Internal array to hold error strings.
char fla_error_string[FLA_MAX_NUM_ERROR_MSGS][FLA_MAX_ERROR_MSG_LENGTH];


/* ***************************************************************************

   FLA_Error_string_for_code()

 *************************************************************************** */

char* FLA_Error_string_for_code( int code )
{
	return fla_error_string[-code];
}

/* ***************************************************************************

   FLA_Error_messages_init()

 *************************************************************************** */

void FLA_Error_messages_init( void )
{
	sprintf( FLA_Error_string_for_code(FLA_INVALID_SIDE),
             "Invalid side parameter value." );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_UPLO),
             "Invalid uplo parameter value." );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_TRANS),
             "Invalid trans parameter value." );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_TRANS_GIVEN_DATATYPE),
	         "Invalid trans value (FLA_CONJ_TRANSPOSE|FLA_CONJ_NO_TRANSPOSE) for given non-complex object datatype" );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_CONJ),
             "Invalid conjugate parameter value." );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_DIRECT),
             "Invalid direction parameter value." );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_STOREV),
             "Invalid storev parameter value." );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_DATATYPE),
             "Invalid datatype value." );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_INTEGER_DATATYPE),
             "Invalid integer datatype value." );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_REAL_DATATYPE),
             "Invalid real datatype value." );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_COMPLEX_DATATYPE),
             "Invalid complex datatype value." );
	sprintf( FLA_Error_string_for_code(FLA_OBJECT_NOT_INTEGER),
             "Expected integer object." );
	sprintf( FLA_Error_string_for_code(FLA_OBJECT_NOT_REAL),
             "Expected real object." );
	sprintf( FLA_Error_string_for_code(FLA_OBJECT_NOT_COMPLEX),
             "Expected complex object." );
	sprintf( FLA_Error_string_for_code(FLA_OBJECT_NOT_SQUARE),
             "Expected square matrix object." );
	sprintf( FLA_Error_string_for_code(FLA_OBJECT_NOT_SCALAR),
             "Expected scalar object." );
	sprintf( FLA_Error_string_for_code(FLA_OBJECT_NOT_VECTOR),
             "Expected vector object." );
	sprintf( FLA_Error_string_for_code(FLA_INCONSISTENT_DATATYPES),
             "Detected inconsistent object datatypes." );
	sprintf( FLA_Error_string_for_code(FLA_NONCONFORMAL_DIMENSIONS),
             "Detected inconsistent object dimensions." );
	sprintf( FLA_Error_string_for_code(FLA_UNEQUAL_VECTOR_DIMS),
             "Detected vectors of unequal dimensions." );
	sprintf( FLA_Error_string_for_code(FLA_NULL_POINTER),
             "Encountered NULL pointer." );
	sprintf( FLA_Error_string_for_code(FLA_SPECIFIED_OBJ_DIM_MISMATCH),
             "Specified dimensions do not match object dimensions." );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_PIVOT_TYPE),
             "Invalid pivot index type specified." );
	sprintf( FLA_Error_string_for_code(FLA_MALLOC_RETURNED_NULL_POINTER),
             "malloc() returned NULL pointer." );
	sprintf( FLA_Error_string_for_code(FLA_OBJECT_BASE_BUFFER_MISMATCH),
             "Detected a buffer address mismatch between adjacent objects." );
	sprintf( FLA_Error_string_for_code(FLA_OBJECTS_NOT_VERTICALLY_ADJ),
             "Object partitions not vertically adjacent." );
	sprintf( FLA_Error_string_for_code(FLA_OBJECTS_NOT_HORIZONTALLY_ADJ),
             "Object partitions not horizontally adjacent." );
	sprintf( FLA_Error_string_for_code(FLA_ADJACENT_OBJECT_DIM_MISMATCH),
             "Object partitions have mismatched dimensions." );
	sprintf( FLA_Error_string_for_code(FLA_OBJECTS_NOT_VERTICALLY_ALIGNED),
             "Object partitions not vertically aligned." );
	sprintf( FLA_Error_string_for_code(FLA_OBJECTS_NOT_HORIZONTALLY_ALIGNED),
             "Object partitions not horizontally aligned." );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_FLOATING_DATATYPE),
             "Expected single or double-precision real or complex datatype value." );
	sprintf( FLA_Error_string_for_code(FLA_OBJECT_NOT_FLOATING_POINT),
             "Expected single or double-precision real or complex object." );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_BLOCKSIZE_VALUE),
             "Detected blocksize value less than or equal to zero." );
	sprintf( FLA_Error_string_for_code(FLA_OPEN_RETURNED_ERROR),
             "open() returned bad file descriptor." );
	sprintf( FLA_Error_string_for_code(FLA_LSEEK_RETURNED_ERROR),
             "lseek() returned error." );
	sprintf( FLA_Error_string_for_code(FLA_CLOSE_RETURNED_ERROR),
             "close() returned error." );
	sprintf( FLA_Error_string_for_code(FLA_UNLINK_RETURNED_ERROR),
             "unlink() returned error." );
	sprintf( FLA_Error_string_for_code(FLA_READ_RETURNED_ERROR),
             "read() returned error." );
	sprintf( FLA_Error_string_for_code(FLA_WRITE_RETURNED_ERROR),
             "write() returned error." );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_QUADRANT),
             "Invalid quadrant parameter value." );
	sprintf( FLA_Error_string_for_code(FLA_NOT_YET_IMPLEMENTED),
             "Function or conditional branch/case not yet implemented." );
	sprintf( FLA_Error_string_for_code(FLA_EXPECTED_NONNEGATIVE_VALUE),
             "Expected non-negative value." );
	sprintf( FLA_Error_string_for_code(FLA_SUPERMATRIX_NOT_ENABLED),
             "SuperMatrix support must be enabled for this code to execute." );
	sprintf( FLA_Error_string_for_code(FLA_UNDEFINED_ERROR_CODE),
             "Undefined error code passed to FLA_Check_error_code()." );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_DIAG),
             "Invalid diag parameter value." );
	sprintf( FLA_Error_string_for_code(FLA_INCONSISTENT_OBJECT_PRECISION),
             "Inconsistent precisions between objects." );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_BLOCKSIZE_OBJ),
             "Encountered blocksize object containing value less than or equal to zero." );
	sprintf( FLA_Error_string_for_code(FLA_VECTOR_DIM_BELOW_MIN),
             "Detected vector dimension below pre-determined minimum." );
	sprintf( FLA_Error_string_for_code(FLA_PTHREAD_CREATE_RETURNED_ERROR),
             "pthread_create() returned error." );
	sprintf( FLA_Error_string_for_code(FLA_PTHREAD_JOIN_RETURNED_ERROR),
             "pthread_join() returned error." );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_ISGN_VALUE),
             "Invalid value for isgn parameter (ie: |isgn| != 1)." );
	sprintf( FLA_Error_string_for_code(FLA_CHOL_FAILED_MATRIX_NOT_SPD),
             "FLA_Chol() failed due to negative diagonal element; matrix not SPD." );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_ELEMTYPE),
             "Invalid object element type value." );
	sprintf( FLA_Error_string_for_code(FLA_POSIX_MEMALIGN_FAILED),
             "posix_memalign() returned error." );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_SUBMATRIX_DIMS),
             "Invalid submatrix dimensions relative to parent matrix." );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_SUBMATRIX_OFFSET),
             "Invalid submatrix offset relative to dimensions of submatrix and parent." );
	sprintf( FLA_Error_string_for_code(FLA_OBJECT_NOT_SCALAR_ELEMTYPE),
             "Object element type is not FLA_SCALAR as expected." );
	sprintf( FLA_Error_string_for_code(FLA_OBJECT_NOT_MATRIX_ELEMTYPE),
             "Object element type is not FLA_MATRIX as expected." );
	sprintf( FLA_Error_string_for_code(FLA_ENCOUNTERED_NON_POSITIVE_NTHREADS),
             "Encountered non-positive (zero) value for number of threads." );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_CONJ_GIVEN_DATATYPE),
	         "Invalid conj value (FLA_CONJUGATE) for given non-complex object datatype" );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_COMPLEX_TRANS),
	         "Invalid complex trans parameter value" );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_REAL_TRANS),
	         "Invalid real trans parameter value" );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_BLAS_TRANS),
	         "Invalid BLAS-style trans parameter value" );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_NONCONSTANT_DATATYPE),
             "Invalid non-constant datatype value." );
	sprintf( FLA_Error_string_for_code(FLA_OBJECT_NOT_NONCONSTANT),
             "Expected non-constant object." );
	sprintf( FLA_Error_string_for_code(FLA_OBJECT_DATATYPES_NOT_EQUAL),
             "Detected unequal object datatypes." );
	sprintf( FLA_Error_string_for_code(FLA_DIVIDE_BY_ZERO),
             "Encountered request to invert zero scalar object." );
	sprintf( FLA_Error_string_for_code(FLA_OBJECT_ELEMTYPES_NOT_EQUAL),
             "Detected unequal object elemtypes." );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_PIVOT_INDEX_RANGE),
             "Invalid pivot index range." );
	sprintf( FLA_Error_string_for_code(FLA_HOUSEH_PANEL_MATRIX_TOO_SMALL),
             "Block-panel Householder matrix is too small." );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_OBJECT_LENGTH),
             "Expected different object length." );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_OBJECT_WIDTH),
             "Expected different object width." );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_ERROR_CHECKING_LEVEL),
             "Invalid error checking level." );
	sprintf( FLA_Error_string_for_code(FLA_ATTEMPTED_OVER_REPART_2X2),
             "Detected attempt to (2x2) repartition more matrix into A11 than exists in source quadrant." );
	sprintf( FLA_Error_string_for_code(FLA_ATTEMPTED_OVER_REPART_2X1),
             "Detected attempt to (2x1) repartition more matrix into A1 than exists in source partition." );
	sprintf( FLA_Error_string_for_code(FLA_ATTEMPTED_OVER_REPART_1X2),
             "Detected attempt to (1x2) repartition more matrix into A1 than exists in source partition." );
	sprintf( FLA_Error_string_for_code(FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED),
             "External LAPACK wrapper was invoked despite not being enabled at configure-time." );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_ROW_STRIDE),
             "Invalid matrix row stride." );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_COL_STRIDE),
             "Invalid matrix column stride." );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_STRIDE_COMBINATION),
             "Invalid combination of matrix row and column strides." );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_VECTOR_DIM),
             "Detected vector of unexpected length/width." );
	sprintf( FLA_Error_string_for_code(FLA_EXPECTED_ROW_VECTOR),
             "Expected object to be a row vector." );
	sprintf( FLA_Error_string_for_code(FLA_EXPECTED_COL_VECTOR),
             "Expected object to be a column vector." );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_INVERSE),
             "Invalid inverse parameter value." );
	sprintf( FLA_Error_string_for_code(FLA_MALLOC_GPU_RETURNED_NULL_POINTER),
             "Attempt to allocate memory on GPU failed." );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_EVD_TYPE),
             "Invalid eigenvalue/vector type parameter value." );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_SVD_TYPE),
             "Invalid singular vector type parameter value." );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_MACHVAL),
             "Invalid machine parameter value." );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_DIAG_OFFSET),
             "Invalid diagonal offset." );
	sprintf( FLA_Error_string_for_code(FLA_EXPECTED_COL_STORAGE),
             "Expected object to be stored by columns." );
	sprintf( FLA_Error_string_for_code(FLA_EXPECTED_ROW_STORAGE),
             "Expected object to be stored by rows." );
	sprintf( FLA_Error_string_for_code(FLA_LAPAC2FLAME_INVALID_RETURN),
             "Invalid return value from lapack2flame interface." );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_SVD_TYPE_COMBINATION),
             "Invalid svd type parameter combination (both parameters are FLA_SVD_VECTORS_OVERWRITE)." );
	sprintf( FLA_Error_string_for_code(FLA_INVALID_SVD_TYPE_AND_TRANS_COMBINATION),
             "Invalid svd type parameters (FLA_SVD_VECTORS_OVERWRITE) and trans parameters combination." );
	sprintf( FLA_Error_string_for_code(FLA_OBJECT_NOT_COMPARABLE),
             "Expected real or int object." );
}

/* ***************************************************************************

   FLA_Print_message()

 *************************************************************************** */

void FLA_Print_message( char *str, char *file, int line )
{
	fprintf( stderr, "\n" );
	fprintf( stderr, "libflame: %s (line %d):\n", file, line );
	fprintf( stderr, "libflame: %s\n", str );
	fflush( stderr );
}

/* ***************************************************************************

   FLA_Abort()

 *************************************************************************** */

void FLA_Abort( void )
{
	fprintf( stderr, "libflame: Aborting.\n");
	//raise( SIGABRT );
	abort();
}

