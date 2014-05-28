function CC=getconstantcoeffs(N,n)
% Given a chebop2, this function converts the problem to one of the form 
% 
%  sum_i  kron(A_i,B_i)
% 
% and returns the discretisation as a cell array 
% 
%  {{   
% 
%       A_1 , B_1 
%       A_2 , B_2
%        .  ,  .
%        .  ,  .
%        .  ,  .
%       A_k , B_k
%
%                   }}
% 

op = N.op;

if isempty(N.coeffs)
maxorder = 2;  % always maxorder = 2 for now. 
A = zeros(maxorder+1);
for jj=0:maxorder
    for kk=0:maxorder
        const = factorial(jj).*factorial(kk);
        f=chebfun2(@(x,y) (x.^jj.*y.^kk)/const);
        A(jj+1,kk+1) = feval(op(f),0,0);
    end
end

% Remove small round off errors
A(abs(A)<10*eps)=0;


N.coeffs=A;  % assign to chebop2 for next time. 
end


% convert matrix of constant coefficients to a discretization using the
% singular value decomposition. 

[na,~]=size(A);
rk = rank(A); [U S V]=svd(A); U = U(:,1:rk); S = S(1:rk,1:rk); V=V(:,1:rk);
CC = cell(rk,2); 
for jj = 1 : rk
    tmp1 = zeros(n); tmp2=zeros(n);
    for kk = 1 : na
         tmp1 = tmp1+U(kk,jj).*spconvermat(n,maxorder-kk+1)*spdiffmat(n,kk-1);
         tmp2 = tmp2+V(kk,jj).*spconvermat(n,maxorder-kk+1)*spdiffmat(n,kk-1);
    end
    CC{jj,1} =tmp1;
    CC{jj,2} =tmp2;
end


end