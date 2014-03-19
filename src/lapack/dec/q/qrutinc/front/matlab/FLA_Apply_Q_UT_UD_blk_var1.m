%
%
%   Copyright (C) 2014, The University of Texas at Austin
%
%   This file is part of libflame and is available under the 3-Clause
%   BSD license, which can be found in the LICENSE file at the top-level
%   directory, or at http://opensource.org/licenses/BSD-3-Clause
%
%

    [ T1T, ...
      T2B ] = FLA_Part_2x1( T1,  b, 'FLA_TOP' );

    W2T = inv( triu( T1T ) )' * ( C1 + D1' * E );

    C1 = C1 - W2T;
    E  = E  - D1 * W2T;

    %------------------------------------------------------------%

    [ DL, DR ] = FLA_Cont_with_1x3_to_1x2( D0, D1, D2, ...
                                           'FLA_LEFT' );

    [ TL, TR ] = FLA_Cont_with_1x3_to_1x2( T0, T1, T2, ...
                                           'FLA_LEFT' );

    [ CT, ...
      CB ] = FLA_Cont_with_3x1_to_2x1( C0, ...
                                       C1, ...
                                       C2, ...
                                       'FLA_TOP' );

  end

  C_out = [ CT
            CB ];

  E_out = E;


return
