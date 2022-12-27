%
%
%   Copyright (C) 2014, The University of Texas at Austin
%
%   This file is part of libflame and is available under the 3-Clause
%   BSD license, which can be found in the LICENSE file at the top-level
%   directory, or at http://opensource.org/licenses/BSD-3-Clause
%
%

function [ B_out ] = FLA_Trsm_rlt_blk_var3( A, B, nb_alg )

  [ BT, ...
    BB ] = FLA_Part_2x1( B, ...
                         0, 'FLA_TOP' );

  while ( size( BT, 1 ) < size( B, 1 ) )

    b = fla_min( size( BB, 1 ), nb_alg );

    [ B0, ...
      B1, ...
      B2 ] = FLA_Repart_2x1_to_3x1( BT, ...
                                    BB, ...
                                    b, 'FLA_BOTTOM' );

    %------------------------------------------------------------%

    B1 = B1 / tril( A )';

    %------------------------------------------------------------%

    [ BT, ...
      BB ] = FLA_Cont_with_3x1_to_2x1( B0, ...
                                       B1, ...
                                       B2, ...
                                       'FLA_TOP' );

  end

  B_out = [ BT
            BB ];

return


