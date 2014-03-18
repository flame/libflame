  function [ AL, AR ] = FLA_Part_1x2( A,...
                                      nb, side )
%
% function [ AL, AR ] = FLA_Part_1x2( A,...
%                                     nb, side )
% Purpose: Partition matrix A into a left and a right side
% where the side indicated by side has nb columns
%
  n = size( A, 2 );
  [ mside, nside ] = size( side );
%
% Check input parameters
%
  if( ( mside ~= 1 )|( nside < 8 )|( nside > 9 ) )
    error('side must be a string with contents equal to FLA_LEFT or FLA_RIGHT');
  elseif( ( nside == 8 )&( ~strcmp( side(1:8), 'FLA_LEFT' ) ) )
    error('side must be a string with contents equal to FLA_LEFT or FLA_RIGHT');
  elseif( nside == 9 )
    if ( ~strcmp( side(1:9), 'FLA_RIGHT' ) )
      error('side must be a string with contents equal to FLA_LEFT or FLA_RIGHT');
    end
  end
%
% Partitioning...
%
  if( strcmp( side(1:8), 'FLA_LEFT' ) )
    AL = A(:,1:nb); 
    AR = A(:,nb+1:n);
  else
    AL = A(:,1:n-nb); 
    AR = A(:,n-nb+1:n);
  end
%
  return;
%
% End of FLA_Part_1x2
%
