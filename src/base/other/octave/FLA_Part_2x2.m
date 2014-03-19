%
%
%   Copyright (C) 2014, The University of Texas at Austin
%
%   This file is part of libflame and is available under the 3-Clause
%   BSD license, which can be found in the LICENSE file at the top-level
%   directory, or at http://opensource.org/licenses/BSD-3-Clause
%
%
  function [ ATL, ATR,...
             ABL, ABR     ] = FLA_Part_2x2( A,...
                                            mb, nb, quadrant )
%
% function [ ATL, ATR,...
%            ABL, ABR     ] = FLA_Part_2x2( A,...
%                                           mb, nb, quadrant )
% Purpose: Partition matrix A into four quadrants
% where the quadrant indicated by quadrant is mb x nb
%
  [ m, n  ] = size( A );
  [ mquadrant, nquadrant ] = size( quadrant );
%
% Check input parameters
%
  if( ( mquadrant ~= 1 )|( nquadrant ~= 6 ) )
    error('quadrant must be a string with contents equal to FLA_TL, FLA_TR, FLA_BL, or FLA_BR');
  elseif( ( ~strcmp( quadrant(1:6), 'FLA_TL' ) )&...
          ( ~strcmp( quadrant(1:6), 'FLA_TR' ) )&...
          ( ~strcmp( quadrant(1:6), 'FLA_BL' ) )&...
          ( ~strcmp( quadrant(1:6), 'FLA_BR' ) ) )
    error('quadrant must be a string with contents equal to FLA_TL, FLA_TR, FLA_BL, or FLA_BR');
  end
%
% Partitioning...
%
  if( strcmp( quadrant(1:6), 'FLA_TL' ) )
    ATL = A(1:mb,1:nb);   
        ATR = A(1:mb,nb+1:n);
    ABL = A(mb+1:m,1:nb); 
        ABR = A(mb+1:m,nb+1:n);
  elseif( strcmp( quadrant(1:6), 'FLA_TR' ) )
    ATL = A(1:mb,1:n-nb);   
        ATR = A(1:mb,n-nb+1:n);
    ABL = A(mb+1:m,1:n-nb); 
        ABR = A(mb+1:m,n-nb+1:n);
  elseif( strcmp( quadrant(1:6), 'FLA_BL' ) )
    ATL = A(1:m-mb,1:nb);   
        ATR = A(1:m-mb,nb+1:n);
    ABL = A(m-mb+1:m,1:nb); 
        ABR = A(m-mb+1:m,nb+1:n);
  else
    ATL = A(1:m-mb,1:n-nb);   
        ATR = A(1:m-mb,n-nb+1:n);
    ABL = A(m-mb+1:m,1:n-nb); 
        ABR = A(m-mb+1:m,n-nb+1:n);
  end
%
  return;
%
% End of FLA_Part_2x2
%
