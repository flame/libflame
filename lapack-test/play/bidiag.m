function [UU,b,VV] = bidiag(a,mode)
%BIDIAG  Upper Bidiagonalization of matrix A.
%        That is,  A = U*B*V' with U and V orthogonal and B
%        having elements on the main diagonal and first upper diagonal.
%        A and B will have the same singular values.
%        A is rectangular and may be real or complex.
%
%usage:  [U,B,V] = bidiag(A)
%
%   or   [U,B,V] = bidiag(A,'imag')
%
%will calculate a complex B if A is complex, otherwise, B will be real
%
%see also: svd, eig, QR, hess, schur, jordan

%Paul Godfrey
%pgodfrey@intersil.com
%Intersil Corp.
%8-12-2002

if exist('mode','var')
   if mode~='imag'
      mode='real';
   end
else
   mode='real';
end

[r,c]=size(a);
n=min(r,c);
UU=eye(r);
VV=eye(c);

for k=1:n

%zero a col
if k<=r
  x=a(k:end,k);
  v=x;

%what if x is a zero vector or has x(1)=0?
  if x(1)==0,
     x(1)=eps*eps;
  end
  v(1)=v(1)+sign(x(1))*norm(x);

  U=eye(length(v))-2/(v'*v)*(v*v');
  z=zeros(k-1,r-k+1);
  U=[eye(k-1) z ;
       z.'    U];
if mode=='real'
  if ~isreal(v(1))
     phi=atan2(imag(v(1)),real(v(1)));
     U(k,:)=exp(-i*phi)*U(k,:);
  end
end
  UU=U*UU;
  a=U*a;
end

%zero a row
if k<c
  x=a(k,k+1:end);
  v=x.';

%what if x is a zero vector or has x(1)=0?
  if x(1)==0,
     x(1)=eps*eps;
  end
  v(1)=v(1)+sign(x(1))*norm(x);

  V=conj(eye(length(v))-2/(v'*v)*(v*v'));
  z=zeros(k,c-k);
  V=[eye(k) z ;
       z.'  V];
if mode=='real'
  if ~isreal(v(1))
     phi=atan2(imag(v(1)),real(v(1)));
     V(:,k+1)=exp(-i*phi)*V(:,k+1);
  end
end
  VV=VV*V;
  a=a*V;
end

end

UU=UU';
if mode=='real'
   a=real(a);

   q=sign(diag(a));
   if r==1,q=q(1:r,1:r);end
   if c==1,q=q(1:c,1:c);end
   q(q>=0)=1;
   q=diag(q);

   if r<=c
      UU=UU*q;
       a=q*a;
   else
      VV=q*VV;
       a=a*q;
   end

end

b=triu(tril((a),1));

if nargout<=1
   UU=b;
end

return