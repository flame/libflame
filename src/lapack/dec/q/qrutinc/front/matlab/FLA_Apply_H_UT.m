function [ a1tout, ...
           A2out ] = FLA_Apply_H_UT( tau, u2, a1t, ...
                                              A2 )
% Compute
% / a1t \  =  / I - 1/tau / 1  \ ( 1  u2^H ) \ / a1t \ 
% \ A2  /     \           \ u2 /             / \ A2  / 
%

  w = ( a1t + u2' * A2 ) / conj( tau );
  
  a1tout = - w + a1t;
  A2out  = - u2 * w + A2;

return
