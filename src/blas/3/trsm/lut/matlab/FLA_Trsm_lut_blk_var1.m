%
%
%   Copyright (C) 2014, The University of Texas at Austin
%
%   This file is part of libflame and is available under the 3-Clause
%   BSD license, which can be found in the LICENSE file at the top-level
%   directory, or at http://opensource.org/licenses/BSD-3-Clause
%
%

function [ B_out ] = FLA_Trsm_lut_blk_var1( A, B, nb_alg )

  [ ATL, ATR, ...
    ABL, ABR ] = FLA_Part_2x2( A, ...
                               0, 0, 'FLA_TL' );

  [ BT, ...
    BB ] = FLA_Part_2x1( B, ...
                         0, 'FLA_TOP' );

  while ( size( ATL, 1 ) < size( A, 1 ) )

    b = fla_min( size( ABR, 1 ), nb_alg );

    [ A00, A01, A02, ...
      A10, A11, A12, ...
      A20, A21, A22 ] = FLA_Repart_2x2_to_3x3( ATL, ATR, ...
                                               ABL, ABR, ...
                                               b, b, 'FLA_BR' );

    [ B0, ...
      B1, ...
      B2 ] = FLA_Repart_2x1_to_3x1( BT, ...
                                    BB, ...
                                    b, 'FLA_BOTTOM' );

    %------------------------------------------------------------%

    B1 = B1 - A01' * B0;
    B1 = triu( A11 )' \ B1;

    %------------------------------------------------------------%

    [ ATL, ATR, ...
      ABL, ABR ] = FLA_Cont_with_3x3_to_2x2( A00, A01, A02, ...
                                             A10, A11, A12, ...
                                             A20, A21, A22, ...
                                             'FLA_TL' );

    [ BT, ...
      BB ] = FLA_Cont_with_3x1_to_2x1( B0, ...
                                       B1, ...
                                       B2, ...
                                       'FLA_TOP' );

  end

  B_out = [ BT
            BB ];

return


