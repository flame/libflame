%
%
%   Copyright (C) 2014, The University of Texas at Austin
%
%   This file is part of libflame and is available under the 3-Clause
%   BSD license, which can be found in the LICENSE file at the top-level
%   directory, or at http://opensource.org/licenses/BSD-3-Clause
%
%

function [ C_out ] = FLA_Syr2k_ln_blk_var5( A, B, C, nb_alg )

  [ AT, ...
    AB ] = FLA_Part_2x1( A, ...
                         0, 'FLA_BOTTOM' );

  [ BT, ...
    BB ] = FLA_Part_2x1( B, ...
                         0, 'FLA_BOTTOM' );

  [ CTL, CTR, ...
    CBL, CBR ] = FLA_Part_2x2( C, ...
                               0, 0, 'FLA_BR' );

  while ( size( AB, 1 ) < size( A, 1 ) )

    b = fla_min( size( AT, 1 ), nb_alg );

    [ A0, ...
      A1, ...
      A2 ] = FLA_Repart_2x1_to_3x1( AT, ...
                                    AB, ...
                                    b, 'FLA_TOP' );

    [ B0, ...
      B1, ...
      B2 ] = FLA_Repart_2x1_to_3x1( BT, ...
                                    BB, ...
                                    b, 'FLA_TOP' );

    [ C00, C01, C02, ...
      C10, C11, C12, ...
      C20, C21, C22 ] = FLA_Repart_2x2_to_3x3( CTL, CTR, ...
                                               CBL, CBR, ...
                                               b, b, 'FLA_TL' );

    %------------------------------------------------------------%

    C21 = C21 + A2 * B1';
    C21 = C21 + B2 * A1';
    C11 = C11 + A1 * B1' + B1 * A1';

    %------------------------------------------------------------%

    [ AT, ...
      AB ] = FLA_Cont_with_3x1_to_2x1( A0, ...
                                       A1, ...
                                       A2, ...
                                       'FLA_BOTTOM' );

    [ BT, ...
      BB ] = FLA_Cont_with_3x1_to_2x1( B0, ...
                                       B1, ...
                                       B2, ...
                                       'FLA_BOTTOM' );

    [ CTL, CTR, ...
      CBL, CBR ] = FLA_Cont_with_3x3_to_2x2( C00, C01, C02, ...
                                             C10, C11, C12, ...
                                             C20, C21, C22, ...
                                             'FLA_BR' );

  end

  C_out = [ CTL, CTR
            CBL, CBR ];

return

