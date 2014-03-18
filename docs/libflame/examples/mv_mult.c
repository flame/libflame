#include "mpi.h"
#include "PLA.h"

int main( int argc, char *argv[] )
{
   int       
      me,  mb, nb,  m, n,  alignCol, alignRow;
   
   MPI_Comm
      comm = MPI_COMM_WORLD;
   
   PLA_Temp          /* PLA_Temp templ = PLA_TEMP_WORLD;          */
      templ;         /* circumvents creating the template,        */
                     /* but restricts specification of a template */
   PLA_Datatype  
      datatype = PLA_DOUBLE;
   
   PLA_Obj
      A, x, y;
                     /* initialize MPI and PLAPACK2e */
   MPI_Init( &argc, &argv );  
   PLA_Init( );               
                     /* get global dimensions and broadcast to all nodes */
   MPI_Comm_rank( comm, &me );
   if ( me == 0 )
   {
      printf( "Enter template blocking size ( mb x nb ): ");
      scanf ( "%d%d", &mb, &nb );
      printf( "Enter matrix dimensions ( m x n ): ");
      scanf ( "%d%d", &m, &n );
      printf( "Enter alignment to template ( i x j ): ");
      scanf ( "%d%d", &alignCol, &alignRow );
   }
   MPI_Bcast( &mb, 1, MPI_INT, 0, comm );  
   MPI_Bcast( &nb, 1, MPI_INT, 0, comm );
   MPI_Bcast( &m,  1, MPI_INT, 0, comm );  
   MPI_Bcast( &n,  1, MPI_INT, 0, comm );
   MPI_Bcast( &alignCol, 1, MPI_INT, 0, comm );  
   MPI_Bcast( &alignRow, 1, MPI_INT, 0, comm );
                     /* create template with recommended values */
   PLA_Temp_create( comm, 1, 1,
		    mb, PLA_ROW_MAJOR, nb, PLA_COL_MAJOR,
		    &templ );
                     /* create distributed matrix A and vectors x and y */
   PLA_Matrix_create( datatype, m, n,
		      alignCol, alignRow,
		      templ, &A );
   PLA_Mcolvector_create( datatype, n, 1,
			  alignRow,
			  templ, &x );
   PLA_Mcolvector_create( datatype, m, 1,
			  alignRow,
			  templ, &y );
                     /* fill matrix A and vectors x and y */
   fill( A );   fill( x );   fill( y );
                     /* perform matrix-vector multiply, y = A x */
   PLA_Gemv( PLA_NO_TRANSPOSE, PLA_ONE, A, x, PLA_ONE, y );
                     /* free objects and template */
   PLA_Obj_free( &A );   PLA_Obj_free( &x );   PLA_Obj_free( &y );
   PLA_Temp_free( &templ );
                     /* finalize PLAPACK2e and MPI */
   PLA_Finalize( );
   MPI_Finalize( );
   return 0;
}
