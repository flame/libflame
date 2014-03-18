#include "PLA.h"

PLA_Error fill( PLA_Obj A )
{
   PLA_Obj 
      AL,   AR,     A0,  a1,  A2,
      aT1,               a01,
      aB1,               alpha11,
                         a21;
   double 
      value;
   
   int 
      i, j, 
      me = PLA_Temp_comm_all_rank( PLA_Obj_template( A ) );
   
   PLA_Part_1x2( A,    
		 &AL,  &AR, 
		 0, PLA_LEFT );
  
   j = 0;
   while ( PLA_Obj_width( AL ) < PLA_Obj_width( A ) )
   {
      PLA_Repart_1x2_to_1x3( AL, /**/ AR,   &A0, /**/ &a1, &A2,
			     1, PLA_RIGHT );
      /*----------------------------------------------------------------*/
      PLA_Part_2x1( a1,   &aT1, 
		          &aB1,
		    0, PLA_TOP );
      i = 0;
      while ( PLA_Obj_length( aT1 ) < PLA_Obj_length( a1 ) )
      {
	 PLA_Repart_2x1_to_3x1( aT1,        &a01, 
			     /* *** */   /* ******** */
				            &alpha11, 
				aB1,        &a21,  
				1, PLA_BOTTOM );
	 /*------------------------------------------------------------*/
	 value = me + i * 0.01 + j * 0.0001;
	 PLA_Obj_set_element( alpha11, &value ); 
	 /*------------------------------------------------------------*/
	 PLA_Cont_with_3x1_to_2x1( &aT1,        a01, 
				                alpha11, 
				/* **** */   /* ******* */
				   &aB1,        a21,
				   PLA_TOP );
	 i++;
      }
      /*------------------------------------------------------------*/
      PLA_Cont_with_1x3_to_1x2( &AL, /**/ &AR,   A0, a1, /**/ A2,
				PLA_LEFT );
      j++;
   }
   
   return PLA_SUCCESS;
}
