%
%
%   Copyright (C) 2014, The University of Texas at Austin
%
%   This file is part of libflame and is available under the 3-Clause
%   BSD license, which can be found in the LICENSE file at the top-level
%   directory, or at http://opensource.org/licenses/BSD-3-Clause
%
%
function [ X ] = Symm_ll_blk_var10( alpha, A, B, C, nb )
%
% Invariant: [ CL, CR ] = [ hatCL, hatCR  + alpha + A * BR]
%
  [ BL, BR ] = FLA_Part_1x2( B, 0, 'FLA_RIGHT' );

  [ CL, CR ] = FLA_Part_1x2( C, 0, 'FLA_RIGHT' );

  while( size( CR, 2 ) ~= size( C, 2 ) )
     b = min( size( CL, 2 ), nb );

    [ C0, C1, C2 ] = FLA_Repart_1x2_to_1x3( CL, CR,
					    b, 'FLA_LEFT' );

    [ B0, B1, B2 ] = FLA_Repart_1x2_to_1x3( BL, BR,
 				            b, 'FLA_LEFT' );
%* ********************************************************************** *%
    C1 = Symm_ll_unb_var10( alpha, A, B1, C1 );
%* ********************************************************************** *%
    [ CL, CR ] = FLA_Cont_with_1x3_to_1x2( C0, C1, C2, 'FLA_RIGHT' );

    [ BL, BR ] = FLA_Cont_with_1x3_to_1x2( B0, B1, B2, 'FLA_RIGHT' );
  end

  X = CR;
  return;
