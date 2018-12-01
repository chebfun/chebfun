function [ZZ, DD, YY] = fadi(A,B,U,V, varargin)
%factored ADI:
%
% fadi(A, B, U, V, p, q)
% solves AX - XB = U*V' in low rank form using factored ADI with 
% ADI shift parameters provided by vectors p, q.
% OUTPUT: ZZ*DD*YY' \approx =  X.
%
% fadi(A, B, U, V, p, q,tol) is as above, except that compression is
% applied as the solution is computed. 
% The compressed solution matches the uncompressed solution
% to a relative accuracy of tol (wrt to the operator norm). 
%
% fadi(A, B, U, V) If A and B have real eigenvalues that are contained in 
% disjoint intervals, the optimal shift parameters are computed automatically
% and the problem is solved to a relative accuracy of approximately machine 
% precision. 
%
% fadi(A, B, U, V, tol) is as above, except the relative accuracy of the 
% of the solution is specified by tol. 
%
% See getshifts_adi and getshifts_smith for help computing shift parameters. 
%
% References: 
%
% [1] Benner, Peter, Ren-Cang Li, and Ninoslav Truhar. 
% "On the ADI method for Sylvester equations." J. of Comp. and App. Math.
% 233.4 (2009): 1035-1045. 

% code written by Heather Wilber (heatherw3521@gmail.com)
% Jan. 2018
%%

[m,r] = size(U);
[n, ~] = size(V);
compute_shifts = 0; 
 
In = speye(n); 
Im = speye(m); 

%parse input 
if isempty(varargin)
    comp = 0;
    compute_shifts = 1;
    tol = eps; 
elseif numel(varargin)==1
    comp = 1; 
    compute_shifts = 1; 
    tol = eps; 
elseif numel(varargin)==2
    comp = 0; 
    p = varargin{1};
    q = varargin{2};
elseif numel(varargin)==3
    comp = 1; 
    p = varargin{1};
    q = varargin{2};
    tol = varargin{3};
end

% user wants shift parameters computed: 
if compute_shifts == 1
    %find intervals where eigenvalues live: 
    a = eigs(A, 1, 'SM'); 
    b = eigs(A, 1, 'LM'); 
    c = eigs(B, 1, 'SM'); 
    d = eigs(B, 1, 'LM'); 
    % determine if the eigenvalues have a complex part. if so, abandon. 
    if any(abs(imag([a b c d])) > 1e-10)
        error('ADI:fadi:cannot automatically compute shift parameters unless the eigenvalues of A and B in AX - XB = F are contained in real, disjoint intervals.')
    end
    % check if intervals overlap
    I1 = [min(a,b) max(a,b)]; 
    I2 = [min(c,d) max(c,d)]; 
    if (( I1(1) < I2(2) && I1(2) > I2(1)) || ( I2(1) < I1(2) && I2(2) > I1(1)) )
        error('ADI:fadi:cannot automatically compute shift parameters unless the eigenvalues of A and B in AX - XB = F are contained in real, disjoint intervals.')
    end
    [p,q] = getshifts_adi([I1 I2], 'tol', tol); 
end

%check for complex-valued shift parameters
Ns = length(p); 
cp = conj(p); 
cq = conj(q); 

B = B';

%factored ADI 
Z(:, 1:r) = (A-Im*q(1))\(U);
Y(:,1:r) = (B-cp(1)*In)\V; 
ZZ = Z;  
YY = Y; 
DD =  (q(1)-p(1)).*ones(r,1); 
    for i = 1:Ns-1 
        Z = Z + (A - q(i+1)*Im)\((q(i+1)-p(i))*Z);
        Y = Y + (B-cp(i+1)*In)\((cp(i+1)-cq(i))*Y);
        ZZ = [ZZ Z]; 
        YY = [YY Y];
        DD = [DD; (q(i+1)-p(i+1)).*ones(r,1)];
    %interim compression step
        if comp==1
            [ZZ, DD, YY] = compression(ZZ, DD, YY, tol, 1);
            DD = diag(DD);  
        end
    end
    DD = diag(DD); 
end


