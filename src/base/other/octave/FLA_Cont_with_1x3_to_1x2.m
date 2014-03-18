  function [ AL, AR ] = FLA_Cont_with_1x3_to_1x2( A0, A1, A2,...
                                                  side )
%
% function [ AL, AR ] = FLA_Cont_with_1x3_to_1x2( A0, A1, A2,...
%                                                 side )
%
% Purpose: Update the 1 x 2 partitioning of matrix A by moving the
% boundaries so that A1 is added to the side indicated by side
%
  [ m0, n0 ] = size( A0 );
  [ m1, n1 ] = size( A1 );
  [ m2, n2 ] = size( A2 );
  [ mside, nside ] = size( side );
%
% Check input parameters
%
  if( ( m0 ~= m1 )|( m1 ~= m2 ) )
    error('input matrices must have the same number of rows');
  elseif( ( mside ~= 1 )|( nside < 8 )|( nside > 9 ) )
    error('side must be a string with contents equal to FLA_LEFT or FLA_RIGHT');
  elseif( ( nside == 8 )&( ~strcmp( side(1:8), 'FLA_LEFT' ) ) )
    error('side must be a string with contents equal to FLA_LEFT or FLA_RIGHT');
  elseif( nside == 9 )
    if ( ~strcmp( side(1:9), 'FLA_RIGHT' ) )
      error('side must be a string with contents equal to FLA_LEFT or FLA_RIGHT');
    end
  end
%
% Continue with...
%
  if( strcmp( side(1:8), 'FLA_LEFT' ) )
    if( ( m0 == 0 )|( n0+n1 == 0 ) )
      AL = zeros( m0, n0+n1 ); 
    else
      AL = [ A0, A1 ];
    end
    AR = A2;
  else
    AL = A0;
    if( ( m0 == 0 )|( n1+n2 == 0 ) )
      AR = zeros( m0, n1+n2 );
    else
      AR = [ A1, A2 ];
    end
  end
%
  return;
%
% End of FLA_Cont_with_1x3_to_1x2
%
