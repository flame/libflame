function [ A0,...
           A1,...
           A2 ] = FLA_Repart_2x1_to_3x1( AT,...
                                         AB,...
                                         mb, side )
%
% Purpose: Repartition a 2 x 1 partitioning of matrix A into
% a 3 x 1 partitioning where submatrix A1 with mb rows is split from
% the side indicated by side
%
  [ mt, nt ] = size( AT );
  [ mbt, nbt ] = size( AB );
  [ mside, nside ] = size( side );

%
% Check input parameters
%
  if( nt ~= nbt )
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
% Repartitioning...
%
  if( strcmp( side(1:7), 'FLA_TOP' ) )
    A0 = AT( 1:mt-mb, : );
    A1 = AT( mt-mb+1:mt, : );
    A2 = AB;
  else
    A0 = AT;
    A1 = AB( 1:mb, : );
    A2 = AB( mb+1:mbt, : );
  end
%
  return;
%
% End of FLA_Repart_2x1_to_3x1
%
