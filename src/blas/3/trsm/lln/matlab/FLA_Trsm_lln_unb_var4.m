
function [ B_out ] = FLA_Trsm_lln_unb_var4( A, B )

  [ BL, BR ] = FLA_Part_1x2( B, ...
                               0, 'FLA_RIGHT' );

  while ( size( BR, 2 ) < size( B, 2 ) )

    [ B0, b1, B2 ]= FLA_Repart_1x2_to_1x3( BL, BR, ...
                                         1, 'FLA_LEFT' );

    %------------------------------------------------------------%

    b1 = tril( A ) \ b1;

    %------------------------------------------------------------%

    [ BL, BR ] = FLA_Cont_with_1x3_to_1x2( B0, b1, B2, ...
                                           'FLA_RIGHT' );

  end

  B_out = [ BL, BR ];

return


