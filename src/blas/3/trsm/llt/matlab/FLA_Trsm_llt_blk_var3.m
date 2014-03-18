
function [ B_out ] = FLA_Trsm_llt_blk_var3( A, B, nb_alg )

  [ BL, BR ] = FLA_Part_1x2( B, ...
                               0, 'FLA_LEFT' );

  while ( size( BL, 2 ) < size( B, 2 ) )

    b = min( size( BR, 2 ), nb_alg );

    [ B0, B1, B2 ]= FLA_Repart_1x2_to_1x3( BL, BR, ...
                                         b, 'FLA_RIGHT' );

    %------------------------------------------------------------%

    B1 = tril( A )' \ B1;

    %------------------------------------------------------------%

    [ BL, BR ] = FLA_Cont_with_1x3_to_1x2( B0, B1, B2, ...
                                           'FLA_LEFT' );

  end

  B_out = [ BL, BR ];

return


