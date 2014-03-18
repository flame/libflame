
#include "FLA_Eig_gest_il.h"
#include "FLA_Eig_gest_iu.h"
#include "FLA_Eig_gest_nl.h"
#include "FLA_Eig_gest_nu.h"

FLA_Error FLA_Eig_gest_internal( FLA_Inv inv, FLA_Uplo uplo, FLA_Obj A, FLA_Obj Y, FLA_Obj B, fla_eig_gest_t* cntl );
FLA_Error FLA_Eig_gest_il( FLA_Obj A, FLA_Obj Y, FLA_Obj B, fla_eig_gest_t* cntl );
FLA_Error FLA_Eig_gest_iu( FLA_Obj A, FLA_Obj Y, FLA_Obj B, fla_eig_gest_t* cntl );
FLA_Error FLA_Eig_gest_nl( FLA_Obj A, FLA_Obj Y, FLA_Obj B, fla_eig_gest_t* cntl );
FLA_Error FLA_Eig_gest_nu( FLA_Obj A, FLA_Obj Y, FLA_Obj B, fla_eig_gest_t* cntl );

