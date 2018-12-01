function [P, Q] = getshifts_fiadi(I,A,B,U,S,V,varargin)
%
% getshifts_fiadi(I,A,B,U,S,V)
% computes a collection of optimal ADI shift parameters for solving 
% AX - XB = U*S*V' with the fiadi method to a relative accuracy approx = 
% machine precision. 
%
% getshifts_fiadi(I,A,B,U,S,V, tol)
% same as above, but the relative accuracy of the solution will be approx. 
% tol. 
%
% getshifts_fiadi(I,A,B,U,S,V, N)
% If N is vector of values, then N(i,:) shift parameters for each
% AX_i - X_iB = F_i used in fiadi are computed (See below:)
%
% In fiadi, the equation AX-XB = F = U*S*V' is split into the following N equations: 
%
%   AX_i - X_iB = F_i, with X_1 + ...+X_N = X, and F_i = U(:,i)*S(i,i)*V(:,i)'. 
%
% The nonzero entries of P(i,:), Q(i,:) are the optimal ADI shift parameters 
% for solving the above eqn. 
%
% See also: getshifts_adi, getshifts_fismith.
% For more help, see Example_lowrankadi12.m
%
%
% References:
% [1] Townsend, Alex, and Heather Wilber. 
% "On the singular values of matrices with high displacement rank." 
% arXiv preprint arXiv:1712.05864 (2017).

% code written by Heather Wilber (heatherw3521@gmail.com)
% Jan, 2018

%%
% part 1: parse input.
s = diag(S); 

if ~isreal(s)
    s = real(s); 
    warning('ADI:getshifts_FIADI',...
    'Nonzero imaginary part of approximate singular values detected. Setting singular values to real part only.')
end

tolmode = 0; 
if isempty(varargin)
    tol = eps; 
    tolmode = 1; 
elseif numel(varargin) == 1
    if isscalar(varargin{1})
    tol = varargin{1}; % compute to tolerance
    tolmode = 1; 
    else
    tolmode = 0;   %compute a fixed # ADI parameters per i.
    N = varargin{1}; 
    Nmax = max(N);
        if ~(length(s)==length(N))
            error('ADI:getshifts_FIADI:number of ADI iterations unspecified for some blocks')
        end
    end
end

%%
% part 2: set up for computing the FIADI shift parameters. 

a = I(1); 
b = I(2); 
c = I(3); 
d = I(4); 
%check if intervals are overlapping: 
I1 = [min(a,b) max(a,b)]; 
I2 = [min(c,d) max(c,d)]; 
if (( I1(1) < I2(2) && I1(2) > I2(1)) || ( I2(1) < I1(2) && I2(2) > I1(1)) )
    error('ADI:getshifts_fiadi: The intervals containing the eigenvalues of A and B in AX-XB = F must be disjoint.')
end

r = length(s); 

[~, Tinv, gam, cr] = mobiusT([a b c d]); 

% figure out number of steps to approx achieve desired tol:
if tolmode==1 
Nmax = ceil(1/pi^2*log(4/tol)*log(16*cr)); %maximum number of iters to ever apply

%estimates req'd for determining tolerance for each block:  
dist = min(abs(b-c), abs(d-a));  
nrmx = estimate_normX(I,A,B,U*sqrt(S),V*sqrt(S)); 
Cst = (4*r/(tol*dist*nrmx)); 

% number of ADI steps for each rank 1 block:

N = min(Nmax, ceil( log(Cst*s)*log(16*cr)/pi^2)); 
N(N<0) =0 ; %if the Cst*s is below tol,  do 0 ADI steps. 
end


Nblk = length(N); 

P = zeros(Nblk, Nmax); % total # blocks by max # iterations.                                          
Q = P;  

%compute shift parameters
if gam > 1e7 
    K = (2*log(2)+log(gam)) + (-1+2*log(2)+log(gam))/gam^2/4;
    m1 = 1/gam^2;        
else
kp = 1-(1/gam)^2; 
K = ellipke(kp); 
end

for j = 1:Nblk
    Ns = N(j); 
    if (gam > 1e7)
        u = (1/2:Ns-1/2)*K/Ns; 
        dn = sech(u) + .25*m1*(sinh(u).*cosh(u)+u).*tanh(u).*sech(u);
    else
    [~, ~, dn] = ellipj((1/2:Ns-1/2)*K/Ns,kp);
    end
    p1 = -gam*dn; %optimal shift parameters Z(T([a b]), T([c,d])). 
    %solve for zeros and poles of rat function on [a b] and [c d]. 
    p = Tinv(p1); 
    q = Tinv(-p1); 
    P(j,1:Ns) = p;
    Q(j,1:Ns) = q;
end
end


function [T, Tinv, gam, M] = mobiusT(I)
%given I = [a b c d] where [a b] [c d] are two disjoint intervals 
% on the real line, T(I) maps to the four points [-gamma, -1, 1, gamma]. 
% M is the cross-ratio.

a = I(1); 
b = I(2); 
c = I(3); 
d = I(4); 

%parameters
M = abs((c-a)*(d-b)/((c-b)*abs(d-a))); 
gam = -1+2*M+2*sqrt(M^2-M); 
A = -gam*a*(1-gam)+gam*(c-gam*d)+c*gam-gam*d; 
B = -gam*a*(c*gam-d)-a*(c*gam-gam*d)-gam*(c*d-gam*d*c); 
C = a*(1-gam)+gam*(c-d)+c*gam-d; 
D = -gam*a*(c-d)-a*(c-gam*d)+c*d-gam*d*c; 

T = @(z) (A*z+B)./(C*z+D);
Tinv = @(z) (D*z-B)./(-C*z+A);
end

function nmx = estimate_normX(I,A,B,U,V) 
%This function applies a few iterations of fADI to estimate ||X||_2, 
% where AX - XB = U*V'.

%first estimate RHS with a rank 1 function:
U = U(:,1); 
V = V(:,1); 

%now apply a few iterations of fADI
[p,q] = getshifts_adi(I, 2); 
[Z, D, Y] = fadi(A, B, U, V, p, q); 

%estimate norm:
nmx = norm(Z)*max(abs(diag(D)))*norm(Y); 
end

