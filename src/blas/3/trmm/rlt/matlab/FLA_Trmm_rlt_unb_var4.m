
function [ B_out ] = FLA_Trmm_rlt_unb_var4( A, B )

  [ BT, ...
    BB ] = FLA_Part_2x1( B, ...
                         0, 'FLA_BOTTOM' );

  while ( size( BB, 1 ) < size( B, 1 ) )

    [ B0, ...
      b1t, ...
      B2 ] = FLA_Repart_2x1_to_3x1( BT, ...
                                    BB, ...
                                    1, 'FLA_TOP' );

    %------------------------------------------------------------%

    b1t = b1t * tril( A )';

    %------------------------------------------------------------%

    [ BT, ...
      BB ] = FLA_Cont_with_3x1_to_2x1( B0, ...
                                       b1t, ...
                                       B2, ...
                                       'FLA_BOTTOM' );

  end

  B_out = [ BT
            BB ];

return


