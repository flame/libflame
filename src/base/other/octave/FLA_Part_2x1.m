  function [ AT,...
             AB ]   = FLA_Part_2x1( A,...
                                    mb, side )
%
% function [ AT,...
%            AB ]   = FLA_Part_2x1( A,...
%                                   mb, side )
%
% Purpose: Partition matrix A into a top and a bottom side
% where the side indicated by {\tt side} has $ {\tt mb} $ rows
%
  m = size( A, 1 );
  [ mside, nside ] = size( side );
%
% Check input parameters
%
  if( ( mside ~= 1 )|( ( nside ~= 7 )&( nside ~= 10 ) ) )
    error('side must be a string with contents equal to FLA_TOP or FLA_BOTTOM');
  elseif( ( nside == 7 )&( ~strcmp( side(1:7), 'FLA_TOP' ) ) )
    error('side must be a string with contents equal to FLA_TOP or FLA_BOTTOM');
  elseif( nside == 10 )
    if ( ~strcmp( side(1:10), 'FLA_BOTTOM' ) )
      error('side must be a string with contents equal to FLA_TOP or FLA_BOTTOM');
    end
  end
%
% Partitioning...
%
  if( strcmp( side(1:7), 'FLA_TOP' ) )
    AT = A(1:mb,:); 
    AB = A(mb+1:m,:);
  else
    AT = A(1:m-mb,:); 
    AB = A(m-mb+1:m,:);
  end
%
  return;
%
% End of FLA_Part_2x1
%
