/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLA_Tevd_n.h"
#include "FLA_Tevd_v.h"

// --- MAC_Tevd_eigval_converged() ---------------------------------------------

#define MAC_Tevd_eigval_converged_ops( eps, safmin, d1, e1, d2 ) \
	fabsf( e1 ) <= (eps) * sqrt( fabsf( d1 ) ) * sqrt( fabsf( d2 ) ) + (safmin)

#define MAC_Tevd_eigval_converged_opd( eps, safmin, d1, e1, d2 ) \
	fabs( e1 )  <= (eps) * sqrt( fabs( d1 ) )  * sqrt( fabs( d2 ) )  + (safmin)

// --- MAC_Tevd_eigval_converged2() ---------------------------------------------

#define MAC_Tevd_eigval_converged2_ops( eps2, safmin, d1, e1, d2 ) \
	(e1) * (e1) <=        (eps2) * fabsf( (d1) * (d2) ) + (safmin)

#define MAC_Tevd_eigval_converged2_opd( eps2, safmin, d1, e1, d2 ) \
	(e1) * (e1) <=        (eps2) * fabs( (d1) * (d2) ) + (safmin)

FLA_Error FLA_Tevd( FLA_Evd_type jobz, FLA_Obj U, FLA_Obj d, FLA_Obj e, FLA_Obj l );

