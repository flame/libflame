
function [ B_out ] = FLA_Trmm_rln_blk_var3( A, B, nb_alg )

  [ BT, ...
    BB ] = FLA_Part_2x1( B, ...
                         0, 'FLA_TOP' );

  while ( size( BT, 1 ) < size( B, 1 ) )

    b = min( size( BB, 1 ), nb_alg );

    [ B0, ...
      B1, ...
      B2 ] = FLA_Repart_2x1_to_3x1( BT, ...
                                    BB, ...
                                    b, 'FLA_BOTTOM' );

    %------------------------------------------------------------%

    B1 = B1 * tril( A );

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


