
function [ C_out ] = FLA_Gemm_nn_blk_var6( alpha, A, B, C, nb_alg )

  [ AL, AR ] = FLA_Part_1x2( A, ...
                             0, 'FLA_RIGHT' );

  [ BT, ...
    BB ] = FLA_Part_2x1( B, ...
                         0, 'FLA_BOTTOM' );

  while ( size( AR, 2 ) < size( A, 2 ) )

    b = min( size( AL, 2 ), nb_alg );

    [ A0, A1, A2 ]= FLA_Repart_1x2_to_1x3( AL, AR, ...
                                         b, 'FLA_LEFT' );

    [ B0, ...
      B1, ...
      B2 ] = FLA_Repart_2x1_to_3x1( BT, ...
                                    BB, ...
                                    b, 'FLA_TOP' );

    %------------------------------------------------------------%

    C = alpha * A1 * B1 + C;

    %------------------------------------------------------------%

    [ AL, AR ] = FLA_Cont_with_1x3_to_1x2( A0, A1, A2, ...
                                           'FLA_RIGHT' );

    [ BT, ...
      BB ] = FLA_Cont_with_3x1_to_2x1( B0, ...
                                       B1, ...
                                       B2, ...
                                       'FLA_BOTTOM' );

  end

  C_out = C;


return

