  function [ AT,...
             AB ] = FLA_Cont_with_3x1_to_2x1( A0,...
                                              A1,...
                                              A2,...
                                              side )
%
% function [ AT,...
%            AB ] = FLA_Cont_with_3x1_to_2x1( A0,...
%                                             A1,...
%                                             A2,...
%                                             side )
%
% Purpose: Update the 2 x 1 partitioning of matrix A by moving the
% boundaries so that A1 is added to the side indicated by side
%
  [ m0, n0 ] = size( A0 );
  [ m1, n1 ] = size( A1 );
  [ m2, n2 ] = size( A2 );
  [ mside, nside ] = size( side );
%
% Check input parameters
%
  if( ( n0 ~= n1 )|( n1 ~= n2 ) )
    error('input matrices must have the same number of columns');
  elseif( ( mside ~= 1 )|( ( nside ~= 7 )&( nside ~= 10 ) ) )
    error('side must be a string with contents equal to FLA_TOP or FLA_BOTTOM');
  elseif( ( nside == 7 )&( ~strcmp( side(1:7), 'FLA_TOP' ) ) )
    error('side must be a string with contents equal to FLA_TOP or FLA_BOTTOM');
  elseif( nside == 10 )
    if ( ~strcmp( side(1:10), 'FLA_BOTTOM' ) )
      error('side must be a string with contents equal to FLA_TOP or FLA_BOTTOM');
    end
  end
%
% Continue with... 
%
  if( strcmp( side(1:7), 'FLA_TOP' ) )
    if( ( m0+m1 == 0 )|( n0 == 0 ) )
      AT = zeros( m0+m1, n0 );
    else
      AT = [ A0;...
             A1 ];
    end
    AB = A2;
  else
    AT = A0;
    if( ( m1+m2 == 0 )|( n0 == 0 ) )
      AB = zeros( m1+m2, n0 );
    else
      AB = [ A1;...
             A2 ];
    end
  end
%
  return;
%
% End of FLA_Cont_with_3x1_to_2x1
%
