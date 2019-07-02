function [p,xi,sumnorms,sumshistory] = watson2(f,n)
%WATSON Best polynomial approximation in the L1-norm for real functions.
% 
% P = WATSON(F, N) computes the best polynomial approximation to the real
% continuous function F in the L1-sense, using Watson's algorithm. F and P 
% are both CHEBFUN objects. 
% 
% Examples:
%   x = chebfun('x'); f = abs(x);
%   p = watson(f, 20); plot(f-p)
% 
% References:
%
%   [1] Watson, G. A. "An algorithm for linear L1 approximation of 
%       continuous functions." IMA Journal of Numerical Analysis, 
%       1.2 (1981): 157-167.
% 
%   [2] Yuji Nakatsukasa and Alex Townsend. "Error localization of best 
%       L1 polynomial approximants." arXiv preprint arXiv:1902.02664 (2019).
%
% See also MINIMAX, CHEBFUN.POLYFIT.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER'S NOTE: Watson's algorithm is a damped Newton iteration that 
% (hopefully) converges to the best L1 polynomial approximation to f. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If f is a polynomial and n > deg(f), then p = f: 
if ( length( f ) <= n+1 &&  numel( domain(f) ) == 2 ) 
    p = f; 
    return
end

% Grab domain of function: 
[a, b] = domain( f ); 

% Compute the polynomial interpolant pn of f at n+1 Chebyshev points. If
% f-pn has roots only at the Chebyshev points, then Theorem 14.4 and 14.5 
% of Powell's Approximation Theory book shows that pn is the best 
% polynomial approximation to f in the L1-norm. 
x = chebpts(n+1, [a, b], 1); 
% [TODO] Compute the Chebyshev interpolant using a DCT for speed:
p_interp = chebfun.interp1(x, feval(f, x), [a b]);  
r = roots( f - p_interp );

if ( numel( r ) == numel( x ) )
    
    % Very easy case, we already have the best L1 polynomial approximant: 
    p = p_interp;
   
else

sums = sign(feval(f,a)-feval(p_interp,a))*intpsign(n,r);
sumshistory = norm(sums);

itwatson = 100; 
for ii = 1:itwatson % Watson update
g = sign(feval(f,a)-feval(p_interp,a))*intpsign(n,r);

T = chebpoly(0:n);
A = feval(T,r);
dfp = diff(f-p_interp);
D = zeros(1,length(r));
for jj = 1:length(r)
    D(jj) = feval(dfp,r(jj));
end
D = diag(2./abs(D)); H = A'*D*A;
% if cond(H)>1e10, keyboard, H = H+norm(H)*eye(size(H))*1e-8; end
dc = H\(g');

dp = chebfun(dc,'coeffs'); % Newton update poly

gam = 1; 
pnow = p_interp; 
while gam>1e-5
p_interp = pnow + gam*dp;
r = roots(f-p_interp);
sums = sign(feval(f,a)-feval(p_interp,a))*intpsign(n,r);
if norm(sums)<sumshistory(end) & length(r)>=n+1, break, end
    gam = gam*.8; 
    % [gam norm(sums)]
end
sumshistory = [sumshistory norm(sums)]; 

if sumshistory(end)<1e-10
    % 'L1min converged'
    break; 
end
end
p = p_interp;
return
end
end

function sums = intpsign(n,f,p)
% find integral of sign(f-p)*T_i
% intpsign(n,f,p) or 
% intpsign(n,r): r are the roots of f-p

if nargin==2,
r = sort(f,'ascend');
else
r = sort(roots(f-p));
end

if length(n)>1,
N = n(2); n = n(1);
else
N = 2*n;
end
[xx,ww] = legpts(N);
sums = ones(1,n);
for ii = 1:n+1
    a = -1; b = r(1); % first interval
    xnow = xx*(b-a)/2+(a+b)/2;
    sums(ii) = ww*cos((ii-1)*acos(xnow))*(b-a);
for jj = 1:length(r)-1
    a = r(jj); b = r(jj+1);
    xnow = xx*(b-a)/2+(a+b)/2;
    sums(ii) = sums(ii) + ((-1)^jj)*ww*cos((ii-1)*acos(xnow))*(b-a);
end
    if length(r)>1
    jj= jj+1;
    a = r(jj); b = 1;
    xnow = xx*(b-a)/2+(a+b)/2;
    sums(ii) = sums(ii) + ((-1)^jj)*ww*cos((ii-1)*acos(xnow))*(b-a);    
    end
end
sums = sums/2; % for correcting the (b-a)/2 term
end

