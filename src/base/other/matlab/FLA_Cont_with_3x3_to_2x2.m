  function [ ATL, ATR,...
             ABL, ABR ] = FLA_Cont_with_3x3_to_2x2( A00, A01, A02,...
                                                    A10, A11, A12,...
                                                    A20, A21, A22, ...
                                                    quadrant )
%
% function [ ATL, ATR,...
%            ABL, ABR ] = FLA_Cont_with_3x3_to_2x2( A00, A01, A02,...
%                                                   A10, A11, A12,...
%                                                   A20, A21, A22, ...
%                                                   quadrant )
%
% Purpose: Update the 2 x 2 partitioning of matrix A by
% moving the boundaries so that A11 is added to
% the quadrant indicated by quadrant
%
  [ m00, n00 ] = size( A00 ); 
      [ m01, n01 ] = size( A01 ); 
          [ m02, n02 ] = size( A02 );
  [ m10, n10 ] = size( A10 ); 
      [ m11, n11 ] = size( A11 ); 
          [ m12, n12 ] = size( A12 );
  [ m20, n20 ] = size( A20 ); 
      [ m21, n21 ] = size( A21 ); 
          [ m22, n22 ] = size( A22 );
  [ mquadrant, nquadrant ] = size( quadrant );
%
% Check input parameters
%
  if( ( m00 ~= m01 )|( m01 ~= m02 ) )
    error('input matrices in first row must have the same number of rows');
  elseif( ( m10 ~= m11 )|( m11 ~= m12 ) )
    error('input matrices in second row must have the same number of rows');
  elseif( ( m20 ~= m21 )|( m21 ~= m22 ) )
    error('input matrices in third row must have the same number of rows');
  elseif( ( n00 ~= n10 )|( n10 ~= n20 ) )
    error('input matrices in first column must have the same number of columns');
  elseif( ( n01 ~= n11 )|( n11 ~= n21 ) )
    error('input matrices in second column must have the same number of columns');
  elseif( ( n02 ~= n12 )|( n12 ~= n22 ) )
    error('input matrices in third column must have the same number of columns');
  elseif( ( mquadrant ~= 1 )|( nquadrant ~= 6 ) )
    error('quadrant must be a string with contents equal to FLA_TL, FLA_TR, FLA_BL, or FLA_BR');
  elseif( ( ~strcmp( quadrant(1:6), 'FLA_TL' ) )&...
          ( ~strcmp( quadrant(1:6), 'FLA_TR' ) )&...
          ( ~strcmp( quadrant(1:6), 'FLA_BL' ) )&...
          ( ~strcmp( quadrant(1:6), 'FLA_BR' ) ) )
    error('quadrant must be a string with contents equal to FLA_TL, FLA_TR, FLA_BL, or FLA_BR');
  end
%
% Continue with...
%
  if( strcmp( quadrant(1:6), 'FLA_TL' ) )
    if( ( m00+m10 == 0 )|( n00+n01 == 0 ) )
      ATL = zeros( m00+m10, n00+n01 );
    else
      ATL = [ A00, A01;...
              A10 A11 ];
    end
    if( ( m02+m12 == 0 )|( n02 == 0 ) )
      ATR = zeros( m02+m12, n02 );
    else
      ATR = [ A02;...
              A12 ];
    end
    if( ( m20 == 0 )|( n20+n21 == 0 ) )
      ABL = zeros( m20, n20+n21 );
    else
      ABL = [ A20, A21 ];
    end
    ABR = A22;
  elseif( strcmp( quadrant(1:6), 'FLA_TR' ) )
    if( ( m00+m10 == 0 )|( n00 == 0 ) )
      ATL = zeros( m00+m10, n00 );
    else
      ATL = [ A00;...
              A10 ];
    end
    if( ( m01+m11 == 0 )|( n01+n02 == 0 ) )
      ATR = zeros( m01+m11, n01+n02 );
    else
      ATR = [ A01, A02;...
              A11, A12 ];
    end
    ABL = A20;
    if( ( m21 == 0 )|( n21+n22 == 0 ) )
      ABR = zeros( m21, n21+n22 );
    else
      ABR = [ A21, A22 ];
    end
  elseif( strcmp( quadrant(1:6), 'FLA_BL' ) )
    if( ( m00 == 0 )|( n00+n01 == 0 ) )
      ATL = zeros( m00, n00+n01 );
    else
      ATL = [ A00, A01 ];
    end
    ATR = A02;
    if( ( m10+m20 == 0 )|( n10+n11 == 0 ) )
      ABL = zeros( m10+m20, n10+n11 );
    else
      ABL = [ A10, A11;...
              A20, A21 ];
    end
    if( ( m12+m22 == 0 )|( n12 == 0 ) )
      ABR = zeros( m12+m22, n12 );
    else
      ABR = [ A12;...
              A22 ];
    end
  else
    ATL = A00;
    if( ( m01 == 0 )|( n01+n02 == 0 ) )
      ATR = zeros( m01, n01+n02 );
    else
      ATR = [ A01, A02 ];
    end
    if( ( m10+m20 == 0 )|( n10 == 0 ) )
      ABL = zeros( m10+m20, n10 );
    else
      ABL = [ A10;...
              A20 ];
    end
    if( ( m11+m21 == 0 )|( n11+n12 == 0 ) )
      ABR = zeros( m11+m21, n11+n12 );
    else
      ABR = [ A11, A12;...
              A21, A22 ];
    end
  end
%
  return;
%
% End of FLA_Cont_with_3x3_to_2x2
%
