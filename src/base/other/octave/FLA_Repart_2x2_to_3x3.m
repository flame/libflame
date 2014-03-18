  function [ A00, A01, A02,...
             A10, A11, A12,...
             A20, A21, A22     ] = FLA_Repart_2x2_to_3x3( ATL, ATR,...
                                                          ABL, ABR,...
                                                          bm, bn, quadrant )
%
% function [ A00, A01, A02,...
%            A10, A11, A12,...
%            A20, A21, A22     ] = FLA_Repart_2x2_to_3x3( ATL, ATR,...
%                                                         ABL, ABR,...
%                                                         bm, bn, quadrant )
%
% Purpose: Repartition a 2 x 2 partitioning of matrix A into
% a 3 x 3 partitioning where mb x nb submatrix A11 is split from
% the quadrant indicated by quadrant
%
  [ mtl, ntl ] = size( ATL ); 
  [ mtr, ntr ] = size( ATR ); 
  [ mbl, nbl ] = size( ABL ); 
  [ mbr, nbr ] = size( ABR ); 
%
% Check input parameters
%
  if( mtl ~= mtr )
    error('input matrices in top row must have the same number of rows');
  elseif( mbl ~= mbr )
    error('input matrices in bottom row must have the same number of rows');
  elseif( ntl ~= nbl )
    error('input matrices in left column must have the same number of columns');
  elseif( ntr ~= nbr )
    error('input matrices in right column must have the same number of columns');
  elseif( ( ~strcmp( quadrant(1:6), 'FLA_TL' ) )&...
         ( ~strcmp( quadrant(1:6), 'FLA_TR' ) )&...
         ( ~strcmp( quadrant(1:6), 'FLA_BL' ) )&...
         ( ~strcmp( quadrant(1:6), 'FLA_BR' ) ) )
    error('quadrant must be a string with contents equal to FLA_TL, FLA_TR, FLA_BL, or FLA_BR');
  end
%
% Repartitioning...
%
  if( strcmp( quadrant(1:6), 'FLA_TL' ) )
    A00 = ATL( 1:mtl-bm, 1:ntl-bn );
        A01 = ATL( 1:mtl-bm, ntl-bn+1:ntl );
            A02 = ATR( 1:mtl-bm, : );
    A10 = ATL( mtl-bm+1:mtl, 1:ntl-bn );
        A11 = ATL( mtl-bm+1:mtl, ntl-bn+1:ntl );
            A12 = ATR( mtl-bm+1:mtl, : );
    A20 = ABL( :, 1:ntl-bn );
        A21 = ABL( :, ntl-bn+1:ntl );
            A22 = ABR;
  elseif( strcmp( quadrant(1:6), 'FLA_TR' ) )
    A00 = ATL( 1:mtr-bm, : );
        A01 = ATR( 1:mtr-bm, 1:bn );
            A02 = ATR( 1:mtr-bm, bn+1:ntr );
    A10 = ATL( mtr-bm+1:mtr, : );
        A11 = ATR( mtr-bm+1:mtr, 1:bn );
            A12 = ATR( mtr-bm+1:mtr, bn+1:ntr );
    A20 = ABL;
        A21 = ABR( :, 1:bn );
            A22 = ABR( :, bn+1:ntr );
  elseif( strcmp( quadrant(1:6), 'FLA_BL' ) )
     A00 = ATL( :, 1:nbl-bn );
         A01 = ATL( :, nbl-bn+1:nbl );
             A02 = ATR;
     A10 = ABL( 1:bm, 1:nbl-bn );
         A11 = ABL( 1:bm, nbl-bn+1:nbl );
             A12 = ABR( 1:bm, : );
     A20 = ABL( bm+1:mbl, 1:nbl-bn );
         A21 = ABL( bm+1:mbl, nbl-bn+1:nbl );
             A22 = ABR( bm+1:mbl, : );
  else
    A00 = ATL;
        A01 = ATR( :, 1:bn );
            A02 = ATR( :, bn+1:nbr );
    A10 = ABL( 1:bm, : );
        A11 = ABR( 1:bm, 1:bn );
            A12 = ABR( 1:bm, bn+1:nbr );
    A20 = ABL( bm+1:mbr, : );
        A21 = ABR( bm+1:mbr, 1:bn );
            A22 = ABR( bm+1:mbr, bn+1:nbr );
  end
%
  return;
%
% End of FLA_Repart_2x2_to_3x3
%
