function [x, w, v, t] = ultrapts(n, lambda, int, meth)
%ULTRAPTS    Ultraspherical points and Gauss-Gegenbauer quadrature weights.
%   [X]= ULTRAPTS(N, LAMBDA) returns N ultraspherical points X in (-1,1)
%   where the ultraspherical weight function is w(x)=(1-x^2)^(LAMBDA-1/2),
%   LAMBDA > -.5, LAMBDA ~=0.
%
%   [X, W] = ULTRAPTS(N, LAMBDA) returns also a row vector W of weights for
%   Gauss-Gegenbauer quadrature.
%   
%   [X, W, V] = ULTRAPTS(N, LAMBDA) returns additionally a column vector V of
%   weights in the barycentric formula corresponding to the points X. The 
%   weights are scaled so that max(abs(V)) = 1.
%
%   ULTRAPTS(N, LAMBDA, INTERVAL) scales the nodes and weights for the 
%   finite interval INTERVAL.
%
%   ULTRAPTS(N, LAMBDA, INTERVAL, METHOD) allows the user to select which 
%   method to use.
%    METHOD = 'REC' uses the recurrence relation for the ultraspherical 
%     polynomials and their derivatives to perform Newton-Raphson iteration. 
%     If 0 < LAMBDA < 1, then the convergence is guaranteed.
%     Default for: LAMBDA <= 3 and N < 100; 3 < LAMBDA <= 8 and N < 500; 
%     8 < LAMBDA <= 13 and N < 1000; 13 < LAMBDA <= 20 and N < 2000;
%     LAMBDA > 20 and N < 3000.
%    METHOD = 'ASY' uses an algorithm adapted from the Hale-Townsend fast 
%     algorithm based upon asymptotic formulae, which is fast and accurate.
%     If 0 < LAMBDA < 1, then the convergence is guaranteed.
%     Default for: LAMBDA <= 3 and N >= 100; 3 < LAMBDA <= 8 and N >= 500; 
%     8 < LAMBDA <= 13 and N >= 1000; 13 < LAMBDA <= 20 and N >= 2000;
%     LAMBDA > 20 and N >= 3000.
%    METHOD = 'GW' will use the traditional Golub-Welsch eigensystem method,
%       which is maintained mostly for historical reasons.
%   
%   [X, W, V, T] = ULTRAPTS(N,LAMBDA) returns also the arccos of the nodes,
%   T = acos(X). In some situations (in particular with 'ASY') these can be
%   computed to a much better relative precision than X.
%
%   The cases LAMBDA=0 and LAMBDA=1 correspond to Gauss-Chebyshev quadratures 
%   nodes and weights, and are treated specially (as a closed form of the nodes 
%   and weights is available). The case LAMBDA=1/2 correspond to Gauss-Legendre
%   quadrature, and it calls LEGPTS which is a more efficient code.
%
% See also CHEBPTS, LEGPTS, JACPTS, LOBPTS, RADAUPTS, HERMPTS, LAGPTS, and
% TRIGPTS.
 
% Copyright 2018 by The Chebfun Developers. See http://www.chebfun.org/ for
% Chebfun information.
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  'GW' by Nick Trefethen, March 2009 - algorithm adapted from [1].
% 'REC' by Lourenco Peixoto, Jul 2015 - see [3].
% 'ASY' by Lourenco Peixoto, Jan 2016 - see [2] and [4].
%
%  References:
%   [1] G. H. Golub and J. A. Welsch, "Calculation of Gauss quadrature rules",
%       Math. Comp. 23:221-230, 1969.
%   [2] N. Hale and A. Townsend, "Fast and accurate computation of Gauss-Legendre 
%       and Gauss-Jacobi quadrature nodes and weights", SIAM J. Sci. Comp., 2013.
%   [3] L. L. Peixoto, "Desigualdades que garantem a convergencia do metodo
%       de Newton-Raphson para os zeros do polinomio ultraesferico no caso
%       principal", Master's thesis, UFMG, Belo Horizonte, 2015.
%   [4] L. L. Peixoto, "Fast, accurate and convergent computation of 
%       Gauss-Gegenbauer quadrature nodes and weights.", In preparation, 2019.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Defaults:
interval = [-1, 1];
method = 'default';
method_set = 0;
 
if (lambda <= -0.5)
    error('CHEBFUN:ultrapts:sizeLAMBDA', 'LAMBDA must be greater than -1/2.')
elseif (lambda >= 30)
    warning('CHEBFUN:ultrapts:largeLAMBDA',...
        'LAMBDA >= 30. Results may not be accurate.')
end
 
% Check inputs:
if ( nargin > 2 )
    if ( nargin == 4 )
        % Calling sequence = ULTRAPTS(N, LAMBDA, INTERVAL, METHOD)
        interval = int;
        method = meth;
        method_set = 1;
    elseif ( nargin == 3 )
        if ( ischar(int) )
            % Calling sequence = ULTRAPTS(N, LAMBDA, METHOD)
            method = int;
            method_set = 1;
        else
            % Calling sequence = ULTRAPTS(N, LAMBDA, INTERVAL)
            interval = int;
        end
    end
    validStrings = {'default', 'GW', 'ASY', 'REC'};
    if ( ~any(strcmpi(method, validStrings)) )
        error('CHEBFUN:ultrapts:inputs', ['Unrecognised input string: ', method]);
    end
    if ( any(isinf(interval)) || numel(interval) ~= 2 || interval(1) >= interval(2) )
        error('CHEBFUN:ultrapts:inputs', 'Interval invalid.');
    end
end
 
 
% Deal with trivial cases:
if (n<0)
    error('CHEBFUN:ultrapts:n', 'First input should be positive number.');
elseif (n==0) % Return empty vectors if n==0
    x = [];
    w = [];
    v = [];
    t = [];
    return
elseif (n==1)
    x = 0;
    w = sqrt(pi)*exp(gammaln(lambda+.5)-gammaln(lambda+1));
    v = 1;
    t = 1;
    [x, w] = rescale(x,w,interval,lambda);
    return
elseif (n==2)
    x = [-1; 1]/sqrt(2*(1+lambda));
    w = sqrt(pi)*exp(gammaln(lambda+.5)-gammaln(lambda+1));
    w = [w, w]/2;
    v = [1 ; -1];
    t = acos(x);
    [x, w] = rescale(x,w,interval,lambda);
    return
end
 
% Special cases:
if ( lambda == 0 ) % Gauss-Chebyshev: lambda = 0
    [x, ~, v] = chebpts(n, interval, 1);
    w = repmat(pi/n,1,n);
    [~, w] = rescale(x, w, interval, lambda);
    t = acos(x);
    return
elseif ( lambda == 1 ) % Gauss-Chebyshev2: lambda = 1
    x = chebpts(n+2, 2);     
    x = x(2:n+1);
    w = pi/(n+1)*(1-x.^2)';
    t = acos(x);
    v = (1-x.^2);  
    v(2:2:end) = -v(2:2:end); 
    [x, w] = rescale(x, w, interval, lambda);
    return
elseif ( lambda == .5 ) % Gauss-Legendre: lambda = 1/2
    [x, w, v, t] = legpts(n, method);
    [x, w] = rescale(x, w, interval, lambda);
    return
end
 
% Choose the method:
t = [];
if ( lambda <= 3 )
    if ( (n < 100 && ~method_set) || strcmpi(method, 'rec') )
        [x, w, v] = rec(n, lambda); % REC (Recurrence relation) [3]
    elseif ( strcmpi(method, 'GW') )
        [x, w, v] = gw(n,lambda);   % GW [1]
    else
        [x, w, v, t] = asy(n,lambda); % HT [2] and Peixoto [4]
    end
elseif ( lambda <= 8 ) 
    if ( (n < 500 && ~method_set) || strcmpi(method, 'rec') )
        [x, w, v] = rec(n, lambda); % REC (Recurrence relation) [3]
    elseif ( strcmpi(method, 'GW') )
        [x, w, v] = gw(n,lambda);   % GW [1]
    else
        [x, w, v, t] = asy(n,lambda); % HT [2] and Peixoto [4]
    end
elseif ( lambda <= 13 )
    if ( (n < 1000 && ~method_set) || strcmpi(method, 'rec') )
        [x, w, v] = rec(n, lambda); % REC (Recurrence relation) [3]
    elseif ( strcmpi(method, 'GW') )
        [x, w, v] = gw(n,lambda);   % GW [1]
    else
        [x, w, v, t] = asy(n,lambda); % HT [2] and Peixoto [4]
    end
elseif ( lambda <= 20 )
    if ( (n < 2000 && ~method_set) || strcmpi(method, 'rec') )
        [x, w, v] = rec(n, lambda); % REC (Recurrence relation) [3]
    elseif ( strcmpi(method, 'GW') )
        [x, w, v] = gw(n,lambda);   % GW see [1]
    else
        [x, w, v, t] = asy(n,lambda); % HT [2] and Peixoto [4]
    end
else % lambda > 20
    if ( (n < 3000 && ~method_set) || strcmpi(method, 'rec') )
        [x, w, v] = rec(n, lambda); % REC (Recurrence relation) [3]
    elseif ( strcmpi(method, 'GW') )
        [x, w, v] = gw(n,lambda);   % GW see [1]
    else
        [x, w, v, t] = asy(n,lambda); % HT [2] and Peixoto [4]
    end
end
    
 
% Compute a T is one is asked for:
if ( nargout == 4 && isempty(t) )
    t = acos(x);
end                
 
% Scale the nodes and quadrature weights:
[x,w] = rescale(x,w,interval,lambda);
 
% Normalise the barycentric weights:
v = abs(v);
v(2:2:end) = -v(2:2:end); 
v = v./max(abs(v));
 
end
 
 
function [x,w] = rescale(x,w,interval,lambda)
% Rescale to arbitrary finite interval:
if ( ~all(interval == [-1 1]) )
    dab = diff(interval);
    x = (x+1)/2*dab + interval(1);
    w = (.5*dab)^(2*lambda)*w;
end
 
end
          
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------Routine for GW algorithm--------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [x, w, v] = gw(n, lambda)
i = (1:n-1)';
bb = .5*sqrt(i.*(i+2*lambda-1)./(i+lambda)./(i+lambda-1));
TT = diag(bb,1) + diag(bb,-1); % Jacobi matrix.
[V, x] = eig( TT ); % Eigenvalue decomposition.
x = diag(x); % Jacobi points.
 
% Quadrature weights:
w = V(1,:).^2*(2^(2*lambda)*exp(2*gammaln(lambda+.5)-gammaln(2*lambda+1))); 
v = sqrt(1-x.^2).*abs(V(1,:))'; % Barycentric weights.
end
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------Routine for REC algorithm-------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [x, w, v] = rec(n,lambda)
 
% Constants:
lam2 = lambda + lambda;
%cte = 2*(n+lambda); % Constant for the weights.
m = floor(.5*(n+1)); % Computes only nonnegative x.
 
% Only for nonnegative x.
k = (1:m).';
theta = (k-(1-lambda)*0.5)/(n+lambda)*pi;
x0 = cos(theta);
cos2 = x0.*x0;
% Sharp initial guess (Forster and Petras, 1993):
x = cos( theta + lambda*(1-lambda)/(2*(n+lambda)^2)*(1-(6+lambda*...
    (1-lambda)*(9-2*cos2))./(12*(n+lambda)^2*...
    (1-cos2))).*cot(theta) );
 
if (lambda>0 && lambda<1)
    % Choose initial guess for guaranteed convergence [3]:
    if (x > x0)
        x = x0;
    end
elseif ( n > 21 )
    % Gatteschi's approximation (1979) for initial guesses on the
    % boundary region:
    Lam = lambda*(1-lambda);
    N = sqrt((n+lambda)^2+lambda*(1-lambda)/3);
    bz = besselroots(lambda-.5,10); % Approximates zeros of Bessel.
    x(1:10) = cos(bz/N-Lam/90.*(bz.^3+2*(lambda^2-lambda-.75).*bz)/N^5);
end
 
dx = inf;
counter = 0;
 
r1=zeros(n,1);
r2=r1;
 
% The terms of recurrence relation for the orthogonal polynomials:
for k = 1:n
    r1(k) = 2*(k-1+lambda)/k;
    r2(k) = (k-2+lam2)/k;
end
 
% Loop until convergence:
while ( norm(dx, inf) > eps && counter < 20 )
    % Initialise:
    P = 1;
    P1 = 0;
    counter = counter + 1;
    for k = 1:n
        P2 = P1;
        P1 = P;
        P = r1(k)*x.*P1 - r2(k)*P2; % P(x) is the orthogonal polynomial.
    end
    PP = (-n*x.*P + (n+lam2-1)*P1)./(1-x.^2);
    % Newton step:
    dx = P./PP;
    % Newton update:
    x = x - dx;
end
 
P = 1;
P1 = 0;
 
 % Once more for P1 and PP:
for k = 1:n
    P2 = P1;
    P1 = P;
    P = r1(k)*x.*P1 - r2(k)*P2; % P(x) is the orthogonal polynomial.
end
c2 = 1-x.*x;
PP = (-n*x.*P + (n+lam2-1)*P1)./c2;
 
% Quadrature weights for nonnegative values:
% w = cte./(c2.*PP.*PP); % Usual relation.
w = 1./(c2.*PP.^2);
  
if ( n >= 20 )
    cte = gammaratio(n+1,lam2-1);
    C = 4^(1-lambda)*pi*cte/exp(2*gammaln(lambda+1));
else
    C = 4^(1-lambda)*pi*exp(gammaln(n+lam2)-gammaln(n+1)-2*gammaln(lambda+1));
end
C = C*lambda*lambda;
w = C*w; 
 
% Reflect for negative values:
s = mod(n,2);
 
x = [-x(1:m-s); x(m:-1:1)];
 
% Reflect for quadrature weights:
w = [w(1:m-s) ; w(m:-1:1)]';
 
% Reflect for derivatives:
ders = [PP(1:m-s); PP(m:-1:1)];
 
% Barycentric weights:
v = 1./ders;
 
end
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------Routine for ASY algorithm------------------------ %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [x,w,v,t] = asy(n,lambda)
% ASY computes nodes and weights using asymptotic formulae.
 
if (n<=21) % Use only interior formula:
    nbdy = floor(.5*(n+1));
    [x, w, v, t] = asy1(n, lambda, nbdy);
    return
end
 
% Determine switch between interior and boundary regions:
nbdy = 10; % Typically, the 10 nodes nearest the boundary.
 
% Interior algorithm:
[x, w, v, t] = asy1(n, lambda, nbdy);
 
% Boundary algorithm:
[x2, w2, v2, t2] = asy2(n, lambda, nbdy);
 
% Combine:
bdyidx1 = n-(nbdy-1):n;
bdyidx2 = nbdy:-1:1;
x(bdyidx1) = x2;
w(bdyidx1) = w2;
v(bdyidx1) = v2;
t(bdyidx1) = t2;
 
% Reflect using symmetry:
x(bdyidx2) = -x2;
w(bdyidx2) = w2;
v(bdyidx2) = v2;
t(bdyidx2) = pi-t2;
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              ASY1 (interior)                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [x,w,v,t] = asy1(n,lambda,nbdy)
% Algorithm for computing nodes and weights in the interior region.
% N > 21.
 
k = (floor((n+1)/2):-1:1);
 
t0 = (k-(1-lambda)*0.5)/(n+lambda)*pi;
cos2 = cos(t0).^2;
% Sharp initial guess (Forster and Petras, 1993):
t = t0 + lambda*(1-lambda)/(2*(n+lambda)^2)*(1-(6+lambda*...
    (1-lambda)*(9-2*cos2))./(12*(n+lambda)^2*...
    (1-cos2))).*cot(t0);
if (lambda>0 && lambda<1)
    % Choose initial guess for guaranteed convergence [4]:
    if (t < t0)
        t = t0;
    end
end
 
 
% Locate the boundary node:
mint = t(end-nbdy+1);
idx = max(find(t<mint,1)-1,1);
 
% Initialise:
dt = inf;
% Newton iteration:
while ( norm(dt, inf) > sqrt(eps)/1000)       % <-- Enough, as again below
    [vals, ders] = feval_asy1(n, t, 0);        % Evaluate via asy formulae
    dt = vals./ders;                           % Newton step
    t = t - dt;
    dt = dt(1:idx-1);       % Ignore boundary terms
end
[vals, ders] = feval_asy1(n, t, 0);  % Once more for good ders. 
t = transpose(t - vals./ders);
 
 
% Constant for the weights: 
% ( cte = 4^(-lambda)*pi*gamma(n+2*lambda)*gamma(n+1)/gamma(n+lambda)^2 )
ratio1 = gammaratio(n+lambda,lambda);
ratio2 = gammaratio(n+lambda,1-lambda);
cte = 4^(-lambda)*pi*ratio1*ratio2;

% Compute x, w, and v:
x = cos(t);
w = (cte./ders.^2);
v = (transpose(sin(t))./ders).';
 
if (mod(n,2))
    x(1) = 0; % computed analytically
    t(1) = pi/2; % computed analytically
    x = [-x(end:-1:2); x];
    w = [w(end:-1:2), w];
    v = -[v(end:-1:2); v];
    t = [pi-t(end:-1:2); t];
else
    x = [-x(end:-1:1); x];
    w = [w(end:-1:1), w];
    v = [-v(end:-1:1); v];
    t = [pi-t(end:-1:1); t];
end
 
 
function [vals, ders] = feval_asy1(n, t, flag) 
% Evaluate asymptotic formula (interior) - Szego p. 197 (8.21.14). 
 
% The polynomial is normalised in order to avoid the computation of 
% gamma(lambda) into the weights:
% If P is the ultraspherical polynomial, and P_norm is the normalised, then
% P_norm = ( gamma(n+1)*gamma(lambda)/(2*gamma(n+lambda)) ) * P.
 
% Number of expansion terms:
r = mod(lambda,2);
if r==0 || r==1
    % If lambda is an integer, then M=lambda-1 gives accurate formula! 
    M = lambda-1;
else
    % M = 30 on otherwise. (Obtained experimentally.)
    M = 30;
end
% Coefficients in expansion:
c = cumprod( (1-lambda:M-lambda)./(n+lambda-1:-1:n+lambda-M) );
d = cumprod((1+lambda:M-1+lambda)./(2:M));
d = lambda*[1, d];
c = [1, c.*d];
% Total number of expansion terms is M+1: 
M=M+1;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some often used vectors/matrices:
onesT = ones(1, length(t));
onesM = ones(M, 1);
Mlam = transpose((0:M-1)+lambda);
Mnlam = transpose(n+lambda-(0:M-1));
onesMcotT = onesM*cot(t);
MnlamonesT = Mnlam*onesT;
MlamonesT = Mlam*onesT;
 
 
twoSinT = onesM*(2*sin(t));
denom = cumprod(twoSinT)./(twoSinT).^(1-lambda);
if ( ~flag )
    alpha = MnlamonesT.*(onesM*t) - .5*pi*MlamonesT;
    cosAlpha = cos(alpha);
    sinAlpha = sin(alpha);
else
    % Adapted from JACPTS, and LEGPTS:
    %%%%%%%%%%%%%%%% Taylor expansion of cos(alpha0) %%%%%%%%%%%%%%
    k = numel(t):-1:1;
    % HI-LO expansion, to accurately compute (n+lambda)*t - (k-.5*lambda)*pi
    ta = double(single(t));
    tb = t - ta;
    hi = n*ta;
    lo = n*tb+lambda*t;
    pia = double(single(pi)); 
    pib = pi - pia;
    dh = (hi-(k-.25)*pia) + lo -.5*(lambda-.5)*pia -...
        (k-.25+.5*(lambda-.5))*pib;
 
    % Compute cosAlpha(1,:) using Taylor series:
    tmp = 0; sgn = 1; fact = 1; DH = dh; dh2 = dh.*dh;
    for jj = 0:100
        dc = sgn*DH/fact;
        tmp = tmp + dc;
        sgn = -sgn;
        fact = fact*(2*jj+3)*(2*jj+2);
        DH = DH.*dh2;
        if ( norm(dc,inf) ) < eps/2000, break, end
    end
    tmp(2:2:end) = -tmp(2:2:end);
    [~, loc] = max(abs(tmp));
    tmp = sign(cos((n+lambda)*t(loc)-.5*lambda*pi)*tmp(loc))*tmp;
    cosAlpha(1,:) = tmp;
 
    % Compute sinAlpha(1,:) using Taylor series:
    tmp = 0; sgn = 1; fact = 1; DH = 1; dh2 = dh.*dh;
    for jj = 0:100
        dc = sgn*DH/fact;
        tmp = tmp + dc;
        sgn = -sgn;
        fact = fact*(2*jj+2)*(2*jj+1);
        DH = DH.*dh2;
        if (norm(dc, inf)) < eps/2000, break, end
    end
    tmp(2:2:end) = -tmp(2:2:end);
    [~, loc] = max(abs(tmp));
    tmp = sign(sin((n+lambda)*t(loc)-.5*lambda*pi)*tmp(loc))*tmp;
    sinAlpha(1,:) = tmp;
 
    % Compute cosAlpha(k,:) and sinAlpha(k,:) for k = 2,...,M:
    sint = sin(t);
    cost = cos(t);
    for kk = 2:M
        cosAlpha(kk,:) = sinAlpha(kk-1,:).*cost - cosAlpha(kk-1,:).*sint;
        sinAlpha(kk,:) = -(cosAlpha(kk-1,:).*cost + sinAlpha(kk-1,:).*sint);
    end
end
% Sum up all the terms:
vals = c*(cosAlpha./denom); % P(theta)
numer = MnlamonesT.*sinAlpha + MlamonesT.*cosAlpha.*onesMcotT;
ders = -c*(numer./denom); % (dP/dtheta)
end
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              ASY2 (boundary)                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [x, w, v, t] = asy2(n, lambda, nbdy) 
% Algorithm for computing nodes and weights near the boundary.
% N > 2.
 
k=(1:nbdy).';
 
% Constant:
Lam = lambda*(1-lambda);
 
% Choose initial guess:
if ( lambda > 1 || lambda < 0 )
    % Gatteschi's approximation (1979):
    N = sqrt((n+lambda)^2+lambda*(1-lambda)/3);
    bz = besselroots(lambda-.5,nbdy); % Approximates zeros of Bessel
    t = bz/N-Lam/90.*(bz.^3+2*(lambda^2-lambda-.75).*bz)/N^5;
else
    % This initial guess guarantees convergence for 0<lambda<1 [4]:
    t0 = (k-(1-lambda)*0.5)/(n+lambda)*pi;
    cos2 = cos(t0).^2;
    % Sharp initial guess (Forster and Petras, 1993):
    t = t0 + lambda*(1-lambda)/(2*(n+lambda)^2)*(1-(6+lambda*...
        (1-lambda)*(9-2*cos2))./(12*(n+lambda)^2*...
        (1-cos2))).*cot(t0);
    % Choose initial guess for guaranteed convergence [4]:
    if (t < t0)
        t = t0;
    end
end
 
% Useful constants:
a=lambda-.5;
b=a;
 
dt = inf; j = 0;
 
cut = 1:10;
% Newton iteration: 
while ( norm(dt,inf) > sqrt(eps)/1000 && j < 20) 
    [tB1, A2, tB2, A3, Ak, Cint, Dint, f2, v] = asy2_firstterms(a, b, t, n);
    Ak=Ak(cut); Cint=Cint(cut,cut); Dint=Dint(cut,cut); f2=f2(cut); v=v(cut); 
    [vals, ders] = feval_asy2(n, t, 0); % Evaluate via asymptotic formula.
    dt = vals./ders;                    % Newton update.
    t = t + dt   ;                      % Next iterate.
    j = j + 1;
end
[tB1, A2, tB2, A3, Ak, Cint, Dint, f2, v] = asy2_firstterms(a, b, t, n);
Ak=Ak(cut); Cint=Cint(cut,cut); Dint=Dint(cut,cut); f2=f2(cut); v=v(cut); 
[vals, ders] = feval_asy2(n, t, 1);     % Evaluate via asymptotic formula.
dt = vals./ders;                        % Newton update
t = t + dt;    
    
% flip:
t = t(nbdy:-1:1); 
ders = ders(nbdy:-1:1);
 
% Constant for the weights:
% (cte = (2^(2*lambda+1)*(n+lambda)^(2*lambda-1)*gamma(n+1)/gamma(n+2*lambda)
ratio = gammaratio(n+2*lambda,1-2*lambda);
cte = 2^(2*lambda+1)*(n+lambda)^(2*lambda-1)*ratio;

% Revert to x-space:
x = cos(t);      
w = (cte./ders.^2).';
if ( find(w==0) )
    warning('CHEBFUN:ultrapts:largeNLAMDBA',...
       'Some WEIGHTS near the boundary region become zero due to overflow.');
end % The overflow occurs on the computation of ders.^2  
 
 
% Constant for the barycentric weights:
% C1 = 2^(2*lambda+.5)/sqrt(pi)*(n+lambda)^(lambda-.5)*gamma(n+lambda)/gamma(n+2*lambda)
ratio = gammaratio(n+2*lambda,-lambda);
C1 = 2^(2*lambda+.5)/sqrt(pi)*(n+lambda)^(lambda-.5)*ratio;

% Revert to x-space:
v = sin(t)./(ders/C1); % barycentric weights with dP/dtheta computed such 
% as in interior region.

% The function below is adapted from JACPTS:
function [vals, ders] = feval_asy2(n, t, flag)
% Evaluate asymptotic formula (boundary) - Baratella and Gatteschi (1998).
% It also computes additional terms for asymptotic series as needed. 
 
% The polynomial is normalised in order to avoid the computation of 
% gamma(lambda) on the weights:
% If P is the ultraspherical polynomial, and P_norm is the normalised, then
% P_norm = ( sqrt(2)*(n+lambda)^(lambda-.5)*gamma(2*lambda)*gamma(n+1) / 
% gamma(lambda+.5)*gamma(n+2*lambda) ) * P.
    
% Useful constants:
rho = n + lambda;
rho2 = rho - 1;
A = Lam;
        
% Evaluate the Bessel functions:
Ja = besselj(a, rho*t, 0);
Jb = besselj(a + 1, rho*t, 0);
Jbb = besselj(a + 1, rho2*t, 0);
if ( ~flag )
    Jab = besselj(a, rho2*t, 0);
else
    % In the final step, perform accurate evaluation
    Jab = besselTaylor(-t, rho*t, a);
end
% Evaluate functions for recursive definition of coefficients:
gt = 2*A*(cot(t)-1./t);
gtdx = .5*A*(4./t.^2-csc(t/2).^2-sec(t/2).^2);
tB0 = .25*gt;
A10 = a*A/6;
A1 = gtdx/8 - (1+2*a)/8*gt./t - gt.^2/32 - A10;
 
% Higher terms:
tB1t = tB1(t); 
A2t = A2(t); 
 
% VALS:
vals = Ja + Jb.*tB0/rho + Ja.*A1/rho^2 + Jb.*tB1t/rho^3 + Ja.*A2t/rho^4;
% DERS:
ders = Jab + Jbb.*tB0/rho2 + Jab.*A1/rho2^2 + Jbb.*tB1t/rho2^3 +...
    Jab.*A2t/rho2^4;
 
% Higher terms:
tB2t = tB2(t); A3t = A3(t);
dv = Jb.*tB2t/rho^5 + Ja.*A3t/rho^6;
vals = vals + dv;
dd = Jbb.*tB2t/rho2^5 + Jab.*A3t/rho2^6;
ders = ders + dd;
   
% Scaling:
valstmp = vals;
denom = (sin(t)/2).^lambda;
vals = (sqrt(t)./denom).*valstmp; % P(theta)
 
% Relation for derivative:
C2 = n/(n+a)*(rho/rho2)^a;
ders = (-n*(2*(n+lambda)-1)*cos(t).*valstmp + 2*(n+a)^2*C2*ders)/...
    (2*(n+lambda)-1);
ders = ders.*(sqrt(t)./(denom.*sin(t))); % dP(theta)
 
% Increase terms as needed:
if ( lambda > 1 || lambda < 0 )
    
    % Initialise:
    del = inf;
    deld = del;
    
    k = 3; % Ak = A3 - A3(1)
    alph = lambda-.5; % constant
    
    while ( del > eps && deld > eps && k<=60 )
        % Akp:
        Akp = Dint*Ak;
        Akp = Akp - Akp(1);
        Akp_t = Akp./t; 
        % Extrapolate point at t = 0:
        w = pi/2-t(2:end);
        w(2:2:end) = -w(2:2:end);
        w(end) = .5*w(end); 
        Akp_t(1) = sum(w.*Akp_t(2:end))/sum(w);
 
        % Bk:
        tBk = -.5*Akp - (.5+alph)*(Cint*Akp_t) + .5*Cint*(f2.*Ak);
        Bk = tBk./t;
        % Extrapolate point at t = 0
        Bk(1) = sum(w.*Bk(2:end))/sum(w);
 
        % Ak1:
        K = Cint*(f2.*tBk);
        Ak1 = .5*(Dint*tBk) - (.5+alph)*Bk - .5*K;
        Ak1 = Ak1 - Ak1(1);
        Atemp = Ak1;
 
        tBk = @(theta) bary(theta,tBk,t,v);
        Ak1 = @(theta) bary(theta,Ak1,t,v);
 
        tBkt = tBk(t); Ak1t = Ak1(t);
 
        dv = Jb.*tBkt/rho^(2*k+1) + Ja.*Ak1t/rho^(2*k+2);
        dd = Jbb.*tBkt/rho2^(2*k+1) + Jab.*Ak1t/rho2^(2*k+2);
 
        del = sqrt(t).*dv./denom;
        vals = vals + del; % P(theta)
 
        deld = ((-n*(2*(n+lambda)-1)*cos(t).*dv + 2*(n+a)^2*C2*dd)...
            /(2*(n+lambda)-1)).*(sqrt(t)./(denom.*sin(t)));
        ders = ders + deld; % dP(theta)
 
        k = k + 1;
        Ak = Atemp;
        del = norm(del./vals,inf);
        deld = norm(deld./ders,inf);
    end
    
end
 
end
 
end
 
% This code is from JACPTS:
function Ja = besselTaylor(t, z, a)
% BESSELTAYLOR    Accurate evaluation of Bessel function J_A for asy2. (See [2].)
% BESSELTAYLOR(T, Z, A) evaluates J_A(Z+T) by a Taylor series expansion about Z. 
 
npts = numel(t);
kmax = 30;
H = bsxfun(@power, t, 0:kmax).';
% Compute coeffs in Taylor expansions about z (See NIST 10.6.7)
[nu, JK] = meshgrid(-kmax:kmax, z);
[Bjk] = besselj(a + nu, JK, 0);
nck = abs(pascal(floor(1.25*kmax), 1)); nck(1,:) = []; % nchoosek
AA = [Bjk(:,kmax+1), zeros(npts, kmax)];
fact = 1;
for k = 1:kmax
    sgn = 1;
    for l = 0:k
        AA(:,k+1) = AA(:,k+1) + sgn*nck(k,l+1)*Bjk(:,kmax+2*l-k+1);
        sgn = -sgn;
    end
    fact = k*fact;
    AA(:,k+1) = AA(:,k+1)/2^k/fact;
end
% Evaluate Taylor series:
Ja = zeros(npts, 1);
for k = 1:npts
    Ja(k,1) = AA(k,:)*H(:,k);
end
end
 
 
%This code is adapted from JACPTS:
function [tB1,A2,tB2,A3,Ak,C,D,f,v] = asy2_firstterms(alph, bet, theta, n)
% ASY2_FIRSTTERMS   First terms for boundary asymptotic series.
% Compute the first terms in asy2 boundary formula.
 
% These constants are more useful than alph and bet:
A = (0.25 - alph^2);
B = (0.25 - bet^2);
 
% For now, just work on half of the domain:
 c = max(max(theta), .5);
if ( n < 30 )
    N = ceil(40 - n);
elseif ( n >= 30 && c > pi/2-.5)
    N = 15;
else
    N = 10;
end
Nm1 = N - 1;
 
 
% Scaled 2nd-kind Chebyshev points and barycentric weights:
t = .5*c*( sin(pi*(-Nm1:2:Nm1)/(2*Nm1)).' + 1 );
v = [.5 ; ones(Nm1,1)];
v(2:2:end) = -1;
v(end) = .5*v(end);
 
% The g's:
g = A*(cot(t/2)  -2./t) - B*tan(t/2);
gp = A*(2./t.^2 - .5*csc(t/2).^2) - .5*(.25-bet^2)*sec(t/2).^2;
gpp = A*(-4./t.^3 + .25*sin(t).*csc(t/2).^4) - 4*B*sin(t/2).^4.*csc(t).^3;
g(1) = 0; gp(1) = -A/6-.5*B; gpp(1) = 0;
 
% B0:
B0 = .25*g./t;
B0p = .25*(gp./t-g./t.^2);
B0(1) = .25*(-A/6-.5*B);
B0p(1) = 0;
 
% A1:
A10 = alph*(A+3*B)/24;
A1 = .125*gp - (1+2*alph)/2*B0 - g.^2/32 - A10;
A1p = .125*gpp - (1+2*alph)/2*B0p - gp.*g/16;
A1p_t = A1p./t;
A1p_t(1) = -A/720 - A^2/576 - A*B/96 - B^2/64 - B/48 + alph*(A/720 + B/48);
 
% Make f accurately: (Taylor series approx for small t)
fcos = B./(2*cos(t/2)).^2;
f = -A*(1/12 + t.^2/240+t.^4/6048 + t.^6/172800 + t.^8/5322240 + ...
    691*t.^10/118879488000 + t.^12/5748019200 + ...
    3617*t.^14/711374856192000 + 43867*t.^16/300534953951232000);
idx = t > .5;
ti = t(idx);
f(idx) = A.*(1./ti.^2 - 1./(2*sin(ti/2)).^2);
f = f - fcos;
 
% Integrals for B1: (Note that N isn't large, so we don't need to be fancy).
C = chebcolloc2.cumsummat(N)*(.5*c);
D = chebcolloc2.diffmat(N)*(2/c);
I = (C*A1p_t);
J = (C*(f.*A1));
 
% B1:
tB1 = -.5*A1p - (.5+alph)*I + .5*J;
tB1(1) = 0;
B1 = tB1./t;
B1(1) = A/720 + A^2/576 + A*B/96 + B^2/64 + B/48 + ...
    alph*(A^2/576 + B^2/64 + A*B/96) - alph^2*(A/720 + B/48);
 
% A2:
K = C*(f.*tB1);
A2 = .5*(D*tB1) - (.5+alph)*B1 - .5*K;
A2 = A2 - A2(1);
 
% A2p:
A2p = D*A2;
A2p = A2p - A2p(1);
A2p_t = A2p./t;
% Extrapolate point at t = 0:
w = pi/2-t(2:end);
w(2:2:end) = -w(2:2:end);
w(end) = .5*w(end);
A2p_t(1) = sum(w.*A2p_t(2:end))/sum(w);
 
% B2:
tB2 = -.5*A2p - (.5+alph)*(C*A2p_t) + .5*C*(f.*A2);
B2 = tB2./t;
% Extrapolate point at t = 0:
B2(1) = sum(w.*B2(2:end))/sum(w);
 
% A3:
K = C*(f.*tB2);
A3 = .5*(D*tB2) - (.5+alph)*B2 - .5*K;
A3 = A3 - A3(1);
 
% temporary term:
Ak = A3;
 
% Make function for output:
tB1 = @(theta) bary(theta, tB1, t, v);
A2 = @(theta) bary(theta, A2, t, v);
tB2 = @(theta) bary(theta, tB2, t, v);
A3 = @(theta) bary(theta, A3, t, v);
end
