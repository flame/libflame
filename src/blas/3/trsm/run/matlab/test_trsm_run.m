%
%
%   Copyright (C) 2014, The University of Texas at Austin
%
%   This file is part of libflame and is available under the 3-Clause
%   BSD license, which can be found in the LICENSE file at the top-level
%   directory, or at http://opensource.org/licenses/BSD-3-Clause
%
%
n = 8;
A = rand( n, n ) + n * eye( n, n );
B = rand( n, n );

Cref = B / triu( A );

norm( FLA_Trsm_run_unb_var1( A, B ) - Cref, 1 )
norm( FLA_Trsm_run_unb_var2( A, B ) - Cref, 1 )
norm( FLA_Trsm_run_unb_var3( A, B ) - Cref, 1 )
norm( FLA_Trsm_run_unb_var4( A, B ) - Cref, 1 )

norm( FLA_Trsm_run_blk_var1( A, B, 3 ) - Cref, 1 )
norm( FLA_Trsm_run_blk_var2( A, B, 3 ) - Cref, 1 )
norm( FLA_Trsm_run_blk_var3( A, B, 3 ) - Cref, 1 )
norm( FLA_Trsm_run_blk_var4( A, B, 3 ) - Cref, 1 )
