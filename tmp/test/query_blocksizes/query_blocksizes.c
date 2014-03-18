
#include "FLAME.h"


int main()
{
    fla_blocksize_t* bp_m;
    fla_blocksize_t* bp_k;
    fla_blocksize_t* bp_n;
    fla_blocksize_t* bp_min;

    FLA_Init();

	bp_m   = FLA_Query_blocksizes( FLA_DIMENSION_M );
	bp_k   = FLA_Query_blocksizes( FLA_DIMENSION_K );
	bp_n   = FLA_Query_blocksizes( FLA_DIMENSION_N );
    bp_min = FLA_Query_blocksizes( FLA_DIMENSION_MIN );

	fprintf( stdout, "             m      k      n    min\n" );
	fprintf( stdout, "float    %5d  %5d  %5d  %5d\n", 
	         FLA_Blocksize_extract( FLA_FLOAT, bp_m ), 
	         FLA_Blocksize_extract( FLA_FLOAT, bp_k ), 
	         FLA_Blocksize_extract( FLA_FLOAT, bp_n ), 
	         FLA_Blocksize_extract( FLA_FLOAT, bp_min ) ); 
	fprintf( stdout, "double   %5d  %5d  %5d  %5d\n", 
	         FLA_Blocksize_extract( FLA_DOUBLE, bp_m ), 
	         FLA_Blocksize_extract( FLA_DOUBLE, bp_k ), 
	         FLA_Blocksize_extract( FLA_DOUBLE, bp_n ), 
	         FLA_Blocksize_extract( FLA_DOUBLE, bp_min ) ); 
	fprintf( stdout, "complex  %5d  %5d  %5d  %5d\n", 
	         FLA_Blocksize_extract( FLA_COMPLEX, bp_m ), 
	         FLA_Blocksize_extract( FLA_COMPLEX, bp_k ), 
	         FLA_Blocksize_extract( FLA_COMPLEX, bp_n ), 
	         FLA_Blocksize_extract( FLA_COMPLEX, bp_min ) ); 
	fprintf( stdout, "dcomplex %5d  %5d  %5d  %5d\n", 
	         FLA_Blocksize_extract( FLA_DOUBLE_COMPLEX, bp_m ), 
	         FLA_Blocksize_extract( FLA_DOUBLE_COMPLEX, bp_k ), 
	         FLA_Blocksize_extract( FLA_DOUBLE_COMPLEX, bp_n ), 
	         FLA_Blocksize_extract( FLA_DOUBLE_COMPLEX, bp_min ) ); 

    FLA_Blocksize_free( bp_m );
    FLA_Blocksize_free( bp_k );
    FLA_Blocksize_free( bp_n );
    FLA_Blocksize_free( bp_min );

    FLA_Finalize();

	return 0;
}

