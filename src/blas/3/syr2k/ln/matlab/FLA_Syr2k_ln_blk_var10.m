
function [ C_out ] = FLA_Syr2k_ln_blk_var10( A, B, C, nb_alg )

  [ AL, AR ] = FLA_Part_1x2( A, ...
                               0, 'FLA_RIGHT' );

  [ BL, BR ] = FLA_Part_1x2( B, ...
                               0, 'FLA_RIGHT' );

  while ( size( AR, 2 ) < size( A, 2 ) )

    b = min( size( AL, 2 ), nb_alg );

    [ A0, A1, A2 ]= FLA_Repart_1x2_to_1x3( AL, AR, ...
                                         b, 'FLA_LEFT' );

    [ B0, B1, B2 ]= FLA_Repart_1x2_to_1x3( BL, BR, ...
                                         b, 'FLA_LEFT' );

    %------------------------------------------------------------%

    C = C + A1 * B1' + B1 * A1';

    %------------------------------------------------------------%

    [ AL, AR ] = FLA_Cont_with_1x3_to_1x2( A0, A1, A2, ...
                                           'FLA_RIGHT' );

    [ BL, BR ] = FLA_Cont_with_1x3_to_1x2( B0, B1, B2, ...
                                           'FLA_RIGHT' );

  end

  C_out = C;


return

