function [ X ] = Symm_ll_unb_var10( alpha, A, B, C )
%
% Invariant: [ CL, CR ] = [ hatCL, hatCR  + alpha + A * BR]
%
  [ BL, BR ] = FLA_Part_1x2( B, 0, 'FLA_RIGHT' );

  [ CL, CR ] = FLA_Part_1x2( C, 0, 'FLA_RIGHT' );

  while( size( CR, 2 ) ~= size( C, 2 ) )
    [ C0, c1, C2 ] = FLA_Repart_1x2_to_1x3( CL, CR,
					    1, 'FLA_LEFT' );

    [ B0, b1, B2 ] = FLA_Repart_1x2_to_1x3( BL, BR,
 				            1, 'FLA_LEFT' );
%* ********************************************************************** *%
    c1 = alpha * ( tril( A ) + tril( A, -1 )' ) * b1 + c1;
%* ********************************************************************** *%
    [ CL, CR ] = FLA_Cont_with_1x3_to_1x2( C0, c1, C2, 'FLA_RIGHT' );

    [ BL, BR ] = FLA_Cont_with_1x3_to_1x2( B0, b1, B2, 'FLA_RIGHT' );
  end

  X = CR;
  return;
