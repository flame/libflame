
function [ B_out ] = FLA_Trsm_lln_blk_var4( A, B, nb_alg )

  [ BL, BR ] = FLA_Part_1x2( B, ...
                               0, 'FLA_RIGHT' );

  while ( size( BR, 2 ) < size( B, 2 ) )

    b = min( size( BL, 2 ), nb_alg );

    [ B0, B1, B2 ]= FLA_Repart_1x2_to_1x3( BL, BR, ...
                                         b, 'FLA_LEFT' );

    %------------------------------------------------------------%

    B1 = tril( A ) \ B1;

    %------------------------------------------------------------%

    [ BL, BR ] = FLA_Cont_with_1x3_to_1x2( B0, B1, B2, ...
                                           'FLA_RIGHT' );

  end

  B_out = [ BL, BR ];

return


