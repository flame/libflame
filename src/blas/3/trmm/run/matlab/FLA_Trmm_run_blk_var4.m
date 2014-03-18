
function [ B_out ] = FLA_Trmm_run_blk_var4( A, B, nb_alg )

  [ BT, ...
    BB ] = FLA_Part_2x1( B, ...
                         0, 'FLA_BOTTOM' );

  while ( size( BB, 1 ) < size( B, 1 ) )

    b = min( size( BT, 1 ), nb_alg );

    [ B0, ...
      B1, ...
      B2 ] = FLA_Repart_2x1_to_3x1( BT, ...
                                    BB, ...
                                    b, 'FLA_TOP' );

    %------------------------------------------------------------%

    B1 = B1 * triu( A );

    %------------------------------------------------------------%

    [ BT, ...
      BB ] = FLA_Cont_with_3x1_to_2x1( B0, ...
                                       B1, ...
                                       B2, ...
                                       'FLA_BOTTOM' );

  end

  B_out = [ BT
            BB ];

return


