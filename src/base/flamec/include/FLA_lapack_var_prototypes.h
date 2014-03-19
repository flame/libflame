/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// Factorizations
#include "FLA_Chol.h"
#include "FLA_LU_nopiv.h"
#include "FLA_LU_piv.h"
#include "FLA_LU_incpiv.h"
#include "FLA_QR_UT.h"
#include "FLA_QR_UT_piv.h"
#include "FLA_QR2_UT.h"
#include "FLA_QR_UT_inc.h"
#include "FLA_LQ_UT.h"
#include "FLA_CAQR2_UT.h"
#include "FLA_CAQR_UT_inc.h"


// Other Decompositions
#include "FLA_Hevd.h"
#include "FLA_Tevd.h"
#include "FLA_Svd.h"
#include "FLA_Bsvd.h"

// Inversions
#include "FLA_Trinv.h"
#include "FLA_SPDinv.h"

// Reductions
#include "FLA_Hess_UT.h"
#include "FLA_Tridiag_UT.h"
#include "FLA_Bidiag_UT.h"

// Solves
#include "FLA_Lyap.h"
#include "FLA_Sylv.h"

// Miscellaneous
#include "FLA_Ttmm.h"
#include "FLA_UDdate_UT.h"
#include "FLA_UDdate_UT_inc.h"

// Utility
#include "FLA_Accum_T_UT.h"
#include "FLA_Apply_G.h"
#include "FLA_Apply_H2_UT.h"
#include "FLA_Apply_HUD_UT.h"
#include "FLA_Apply_Q_UT.h"
#include "FLA_Apply_Q2_UT.h"
#include "FLA_Apply_CAQ2_UT.h"
#include "FLA_Apply_QUD_UT.h"
#include "FLA_Apply_Q_UT_inc.h"
#include "FLA_Apply_CAQ_UT_inc.h"
#include "FLA_Apply_QUD_UT_inc.h"
#include "FLA_Apply_pivots.h"

// Eigensolvers
#include "FLA_Eig_gest.h"
