/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// --- MAC_Bsvd_sinval_is_converged() ------------------------------------------

#define MAC_Bsvd_sinval_is_converged_ops( tol, d1, e1 ) \
	fabsf( (e1) ) <= fabsf( (tol) * (d1) )

#define MAC_Bsvd_sinval_is_converged_opd( tol, d1, e1 ) \
	fabs(  (e1) ) <= fabs(  (tol) * (d1) )

// --- FLA_Bsvd_sinval_v_opt_var1() --------------------------------------------

FLA_Error FLA_Bsvd_sinval_v_opt_var1( FLA_Obj tol, FLA_Obj thresh, FLA_Obj G, FLA_Obj H, FLA_Obj d, FLA_Obj e, FLA_Obj n_iter );
FLA_Error FLA_Bsvd_sinval_v_ops_var1( integer       m_A,
                                      integer       n_GH,
                                      integer       n_iter_allowed,
                                      float     tol, 
                                      float     thresh, 
                                      scomplex* buff_G, integer rs_G, integer cs_G,
                                      scomplex* buff_H, integer rs_H, integer cs_H,
                                      float*    buff_d, integer inc_d, 
                                      float*    buff_e, integer inc_e,
                                      integer*      n_iter );
FLA_Error FLA_Bsvd_sinval_v_opd_var1( integer       m_A,
                                      integer       n_GH,
                                      integer       n_iter_allowed,
                                      double    tol, 
                                      double    thresh, 
                                      dcomplex* buff_G, integer rs_G, integer cs_G,
                                      dcomplex* buff_H, integer rs_H, integer cs_H,
                                      double*   buff_d, integer inc_d, 
                                      double*   buff_e, integer inc_e,
                                      integer*      n_iter );

