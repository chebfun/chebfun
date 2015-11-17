function [x, w, v, t] = ultrapts(n, lambda, int, meth, conv)
%ULTRAPTS    Ultraspherical points and Gauss-Gegenbauer quadrature weights.
%   [X]= ULTRAPTS(N, LAMBDA) returns N ultraspherical points X in (-1,1) with
%   parameter 0<=LAMBDA<=1 where the ultraspherical weight function is
%   defined by w(x)=(1-x^2)^(LAMBDA-1/2).
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
%   ULTRAPTS(N, LAMBDA, INTERVAL, METHOD) or ULTRAPTS(N, LAMBDA, METHOD) allows 
%   the user to select which method to use.
%    METHOD = 'REC' uses the recurrence relation for the ultraspherical polynomials
%     and their derivatives to perform Newton-Raphson iteration with some
%     approximations that guarantee convergence to the roots - see [3]. 
%     Default for N < 100.
%    METHOD = 'ASY' uses the Hale-Townsend fast algorithm based upon asymptotic
%     formulae, which is fast and accurate - see [2]. Default for N >= 100. 
%    METHOD = 'GW' will use the traditional Golub-Welsch eigensystem method,
%       which is maintained mostly for historical reasons - see [1].
%   
%   ULTRAPTS(N, LAMBDA, INTERVAL, 'ASY', CONVERGENCE) or 
%   ULTRAPTS(N, LAMBDA, 'ASY', CONVERGENCE) with CONVERGENCE ~=0 uses 'ASY'
%   with some approximations that guarantee convergence to the roots - see
%   [4]. Default is CONVERGENCE=0. The parameter CONVERGENCE with 'REC' or
%   'GW' doesn't have any meaning.
%
%   [X, W, V, T] = ULTRAPTS(N,LAMBDA) returns also the arccos of the nodes,
%   T = acos(X). In some situations (in particular with 'ASY') these can be
%   computed to a much better relative precision than X.
%
%   The cases LAMBDA=0 and LAMBDA=1 correspond to Gauss-Chebyshev quadratures 
%   nodes and weights, and are treated specially (as a closed form of the nodes 
%   and weights is available). The case LAMBDA=.5 correspond to Gauss-Legendre
%   quadrature, and it calls LEGPTS which is a more efficient code.
%
% See also CHEBPTS, LEGPTS, JACPTS, LOBPTS, RADAUPTS, HERMPTS, LAGPTS, and
% TRIGPTS.

% Copyright 2015 by The Chebfun Developers. See http://www.chebfun.org/ for
% Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  'GW' by Nick Trefethen, March 2009 - algorithm adapted from [1].
% 'REC' by Nick Hale, July 2011
% 'ASY' by Nick Hale & Alex Townsend, May 2012 - see [2].
%
%
%  References:
%   [1] G. H. Golub and J. A. Welsch, "Calculation of Gauss quadrature rules",
%       Math. Comp. 23:221-230, 1969,
%   [2] N. Hale and A. Townsend, "Fast and accurate computation of Gauss-Legendre 
%       and Gauss-Jacobi quadrature nodes and weights", SIAM J. Sci. Comp., 2013.
%   [3] L. L. Peixoto, "Desigualdades que garantem a convergência do método
%       de Newton-Raphson para os zeros do polinômio ultraesférico no caso
%       principal", Master's thesis, UFMG, Belo Horizonte, 2015.
%   [4] L. L. Peixoto, "On the convergence of Newton-Raphson to the zeros of
%       ultraspherical polynomials on theta-variable", In preparation, 2016.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Defaults:
interval = [-1, 1];
method = 'ASY';
convergence = 0;
method_set = 0;

if (lambda < 0 || lambda > 1)
    error('CHEBFUN:ultrapts:lambda', '0<=LAMBDA<=1'); 
    % [TODO]: lambda<0 and lambda>1.
elseif (lambda<30*eps && lambda~=0) % ASY is not accurate for 0<LAMBDA<30*eps.
    warning('CHEBFUN:ultrapts:smallLAMBDA',...
        '0<LAMBDA<30*eps.  Results may not be accurate with ASY.')
end % [TODO]: What is the smallest NUMBER<1 for which ASY is not accurate? 
    % NUMBER=.99? 

% Check the inputs:
if ( nargin > 2 )
    if ( nargin == 5 )
        % Calling sequence = ULTRAPTS(N, LAMBDA, INTERVAL, METHOD=ASY, CONV)
        interval = int;
        method_set=1;
        convergence = conv;
    elseif ( nargin == 4 )
        if ( ischar(int) ) % Calling sequence = ULTRAPTS(N, LAMBDA, METHOD=ASY, CONV)
            method = int;
            convergence = meth;
            method_set=1;
        else % Calling sequence = ULTRAPTS(N, LAMBDA, INTERVAL, METHOD)
            interval = int;
            method = meth;
            method_set=1;
        end
    elseif ( nargin == 3 )
        if ( ischar(int) ) % Calling sequence = ULTRAPTS(N, LAMBDA, METHOD)
            method = int;
            method_set=1;
        else % Calling sequence = ULTRAPTS(N, LAMBDA, INTERVAL)
            interval = int;
        end
    end
    validStrings = {'default', 'GW', 'asy', 'rec'};
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
    w = gamma(lambda+.5)*sqrt(pi)/gamma(lambda+1);
    v = 1;
    t = 1;
    [x, w] = rescale(x,w,interval,lambda);
    return
elseif (n==2)
    x = [-1; 1]/sqrt(2*(1+lambda));
    w = gamma(lambda+.5)*sqrt(pi)/gamma(lambda+1);
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
elseif ( lambda == 1)   % Gauss-Chebyshev2: lambda = 1
    x = chebpts(n+2, 2);     
    x = x(2:n+1);
    w = pi/(n+1)*(1-x.^2)';
    t = acos(x);
    v = (1-x.^2);  
    v(2:2:end) = -v(2:2:end); 
    [x, w] = rescale(x,w,interval,lambda);
    return
elseif (lambda == .5) % Gauss-Legendre: lambda = 1/2
    [x, w, v, t] = legpts(n, interval);
    return
end

% Choose the method:
t = [];
if ( (n < 100 && ~method_set) || strcmpi(method, 'rec') )
    [x, w, v] = rec(n,lambda);  % REC with convergence guaranteed see [3]
elseif ( strcmpi(method, 'GW') )
    [x, w, v] = gw(n,lambda);   % GW see [1]
else
    [x, w, v, t] = asy(n,lambda,convergence); % HT see [2]. For convergence
    % see Peixoto [4]
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
w = V(1,:).^2*(2^(2*lambda)*gamma(lambda+.5)^2/gamma(2*lambda+1)); 
v = sqrt(1-x.^2).*abs(V(1,:))'; % Barycentric weights.
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------Routine for REC algorithm-------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, w, v] = rec(n,lambda)

% Constant for w: (cte = 4^(1-lambda)*pi*gamma(n+2*lambda)/gamma(n+1)/gamma(lambda)^2)
cte = 4^(1-lambda)*pi*exp(gammaln(n+2*lambda)-gammaln(n+1)-2*gammaln(lambda));

% Constant:
m = floor(.5*n);

if n > 1
    % Only for positive x.
    k = (1:m).';
    theta = (k-(1-lambda)*.5)/(n+lambda)*pi;
    % This initial guess guarantees convergence [3]:
    x0 = cos(theta);
    cos2 = x0.*x0;
    % Sharp initial guess (Förster and Petras, 1993):
    x = cos(theta + lambda*(1-lambda)/(2*(n+lambda)^2)*(1-(6+lambda*...
        (1-lambda)*(9-2*cos2))./(12*(n+lambda)^2*...
        (1-cos2))).*cot(theta));
    % Choose initial guess that guarantees convergence [3]:
    if  x > x0 
        x = x0;
    end
    
    % Initialise:
    Pm2 = 1;
    Pm1 = 2*lambda*x;
    dx = inf;
    counter = 0;
    
    r1=zeros(n,1);
    r2=r1;

    % The terms of recurrence relation:
    for k = 1:n-1
        r1(k) = 2*(k+lambda);
        r2(k) = (1-k-2*lambda);
    end
        
    % Loop until convergence:
    while ( norm(dx, inf) > eps && counter < 10 )
        counter = counter + 1;
        for k = 1:n-1,
            P = (r1(k)*Pm1.*x+r2(k)*Pm2)/(k+1);
            Pm2 = Pm1;
            Pm1 = P;
        end
        PP = ((n+2*lambda-1)*Pm2-n*x.*Pm1)./(1-x.*x);
        % Newton step:
        dx = P./PP;
        % Newton update:
        x = x - dx;
        % Reinitialise:
        Pm2 = 1;
        Pm1 = 2*lambda*x;
    end

    % Once more for Pm1 and PP (Required for Yakimiw's formula). 
    for k = 1:n-1,
        P = (r1(k)*Pm1.*x+r2(k)*Pm2)/(k+1);
        Pm2 = Pm1;
        Pm1 = P;
    end
    c2 = 1-x.*x;
    PP = ((n+2*lambda-1)*Pm2-n*x.*Pm1)./c2;

    % Quadrature weights for positive values:
    % w = cte./(c2.*PP.*PP); % Usual relation. 
    w = cte*c2./(c2.*PP-2*lambda*x.*Pm1+.5*(n*(n+2*lambda)+2*lambda...
        ./c2).*Pm1.^2./PP).^2; % Yakimiw's formula (more accurate than 
    % usual relation).
end

% If n is odd then computes x=0 and corresponding w analytically:
if ( mod(n,2) )
    s = 1;
    x(m+1,1) = 0;
    % Calculate PP analytically: (PP =
    % (n+2*lambda-1)*gamma((n-1)/2+lambda)/gamma((n+1)/2)/gamma(lambda);
    PP(m+1,1) = 2*exp(gammaln((n+1)/2+lambda)-gammaln((n+1)/2)-gammaln(lambda));
    PP2 = PP(m+1,1)*PP(m+1,1); % PP^2
    w(m+1,1) = cte/PP2;
else
    s = 0;
end

% Reflect for negative values:
x = [-x(1:m) ; x(m+s:-1:1)];

% Reflect for quadrature weights:
w = [w(1:m) ; w(m+s:-1:1)]';

% Reflect for derivatives:
ders = [PP(1:m); PP(m+s:-1:1)];

% Barycentric weights:
v = 1./ders;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------------Routine for ASY algorithm------------------------ %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,w,v,t] = asy(n,lambda,convergence)
% ASY computes nodes and weights using asymptotic formulae. If 
% CONVERGENCE ~=0 then the convergence is guaranteed [4].

% Constant for w:
% ( cte = 4^(1-lambda)*pi*gamma(n+2*lambda)/gamma(n+1)/gamma(lambda)^2 )
% See ratio of gamma functions in [2]. 
ds = .5*(2*lambda-1)^2/n;
s0 = ds;
j = 1; 
while ( abs(ds/s0) > eps/100 ) % Taylor series in expansion 
    j = j+1;
    ds = -(2*lambda-1)*(j-1)/(j+1)/n*ds;
    s0 = s0 + ds;
end
p2 = exp(s0)*sqrt(n+2*lambda-1)*n^(2*lambda-1.5);
% Stirling's series:
g = [1, 1/12, 1/288, -139/51840, -571/2488320, 163879/209018880, ...
    5246819/75246796800, -534703531/902961561600, ...
    -4483131259/86684309913600, 432261921612371/514904800886784000];
f = @(z) sum(g.*[1, cumprod(ones(1, 9)./z)]);
cte = 4^(1-lambda)*pi/gamma(lambda)^2*p2*(f(n+2*lambda-1)/f(n));

if (n<=21) % Use only boundary formula:
    nbdy = floor(.5*(n+1));
    [x2, w2, v2, t2] = asy2(n, lambda, nbdy, convergence, cte);
    s = mod(n,2);
    
    x = [-x2(end:-1:1+s); x2(1:end)];
    w = [w2(end:-1:1+s); w2(1:end)]';
    v = [-v2(end:-1:1+s); v2(1:end)];
    t = [pi-t2(end:-1:1+s); t2(1:end)];
    return
end

% Determine switch between interior and boundary regions:
nbdy = 10; % Typically, the 10 nodes nearest the boundary.
% Interior algorithm:
[x1, w1, v1, t1] = asy1(n, lambda, nbdy, convergence, cte);
% Boundary algorithm:
[x2, w2, v2, t2] = asy2(n, lambda, nbdy, convergence, cte);
% Combine:
s = mod(n,2);
m = floor((n+1)/2);
bdy1 = (10:-1:1);
bdy2 = (1:10);
int1 = (m-10:-1:1+s);
int2 = (1:m-10);
x = [-x2(bdy1); -x1(int1); x1(int2); x2(bdy2)];
w = [w2(bdy1); w1(int1); w1(int2); w2(bdy2)]';
v = [v2(bdy1); -v1(int1); -v1(int2); v2(bdy2)];
t = [pi-t2(bdy1); pi-t1(int1); t1(int2); t2(bdy2)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              ASY1 (interior)                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,w,v,t] = asy1(n,lambda,nbdy,convergence,cte)
% Algorithm for computing nodes and weights in the interior region.
% N > 21.
% CONVERGENCE ~=0 garantees convergence of Newton-Raphson [4].

k = (floor((n+1)/2):-1:1);

if (convergence ~=0) 
    % This initial guess guarantees convergence [4]:
    t0 = (k-(1-lambda)*0.5)/(n+lambda)*pi;
    cos2 = cos(t0).^2;
    % Sharp initial guess (Förster and Petras, 1993):
    t = t0 + lambda*(1-lambda)/(2*(n+lambda)^2)*(1-(6+lambda*...
        (1-lambda)*(9-2*cos2))./(12*(n+lambda)^2*...
        (1-cos2))).*cot(t0);
    % Choose initial guess that guarantees convergence [4]:
    if  t < t0 
        t = t0;
    end
else
    % Gatteschi's approximation (1979):
    LAMBDA = lambda*(1-lambda);
    N = sqrt((n+lambda)^2+LAMBDA);
    u = (n+1-2*k)*pi/(2*N);
    tu = tan(u); tu2 = tu.*tu;
    A1 = LAMBDA*(tu-u);
    A2 = 3*LAMBDA^2*u.*(1-2*tu2)+LAMBDA*tu.*(6-3*LAMBDA+(6+7*LAMBDA).*...
        tu2);
    x = sin(u-.5*A1/N^2+A2/(24*N^4));
    t = acos(x);
end

% Locate the boundary node:

mint = t(end-nbdy+1);
idx = max(find(t<mint,1)-1,1);

% Initialise:
dt = inf;
% Newton iteration: 
while ( norm(dt, inf) > sqrt(eps)/1000 )       % <-- Enough, as again below
    [vals, ders] = feval_asy1(n, t, 0);        % Evaluate via asy formulae
    dt = vals./ders;                           % Newton step
    t = t - dt;
    dt = dt(1:idx-1);                          % Ignore boundary terms
end
[vals, ders] = feval_asy1(n, t, 1);  % Once more for good ders. 

t = transpose(t - vals./ders);

% Compute x, w, and v:
x = cos(t);
w = (cte./ders.^2).';
v = (transpose(sin(t))./ders).';

if (mod(n,2))
    x(1)=0;

    % Calculate ders analytically: (ders =
    % (n+2*lambda-1)*gamma((n-1)/2+lambda)/gamma((n+1)/2)/gamma(lambda)
    % See ratio of gamma functions in [2]. 
    k = .5*(n-1);
    ds = .5*(lambda-1)^2/k;
    s0 = ds;
    j = 1;
    while ( abs(ds/s0) > eps/100 ) % Taylor series in expansion 
        j = j+1;
        ds = -(lambda-1)*(j-1)/(j+1)/k*ds;
        s0 = s0 + ds;
    end
    p2 = exp(s0)*sqrt(k+lambda-1)*k^(lambda-1.5);
    % Stirling's series:
    g = [1, 1/12, 1/288, -139/51840, -571/2488320, 163879/209018880, ...
        5246819/75246796800, -534703531/902961561600, ...
        -4483131259/86684309913600, 432261921612371/514904800886784000];
    f = @(z) sum(g.*[1, cumprod(ones(1, 9)./z)]);
    ders = (n+2*lambda-1)*p2*(f(k+lambda-1)/f(k))/gamma(lambda);
    ders2 = ders*ders; % ders^2
    
    w(1) = cte/ders2;
    v(1) = 1/ders;
    t(1) = pi/2;
end


function [vals, ders] = feval_asy1(n, t, flag) 
% Evaluate asymptotic formula (interior) - Szegö p. 197 (8.21.15).

% Max number of expansion terms:
M = 20; % Assuming similar to LEGPTS.

% Coefficients in expansion:
c = cumprod( (lambda:lambda+M-1)./(1:M) );
d = cumprod( (1-lambda:M-lambda)./(n+lambda+1:n+lambda+M) );
c = [1, c.*d];

% Constant out the front: 
% (C = 2*sin(lambda*pi)/pi*gamma(n+2*lambda)*gamma(1-lambda)/gamma(n+lambda+1))
% See ratio of gamma functions in [2]. 
ds = .5*(lambda-1)^2/(n+lambda);
s = ds;
j = 1;
while ( abs(ds/s) > eps/100 ) % Taylor series in expansion 
    j = j+1;
    ds = -(lambda-1)*(j-1)/(j+1)/(n+lambda)*ds;
    s = s + ds;
end
p2 = exp(s)*sqrt(n+2*lambda-1)*(n+lambda)^(lambda-1.5);
% Stirling's series:
g = [1, 1/12, 1/288, -139/51840, -571/2488320, 163879/209018880, ...
    5246819/75246796800, -534703531/902961561600, ...
    -4483131259/86684309913600, 432261921612371/514904800886784000];
f = @(z) sum(g.*[1, cumprod(ones(1, 9)./z)]);
C = 2*sin(lambda*pi)/pi*gamma(1-lambda)*p2*(f(n+2*lambda-1)/f(n+lambda));

% How many terms required in the expansion? (Avoid constant C because M=1 
% when lambda tends to 0).
vec = (lambda:M+1);
vec = vec(1:21);
R = c.*2./(2*sin(mint)).^(vec);
R = R(abs(R) > eps);
M = length(R);
c = c(1:M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Some often used vectors/matrices:
onesT = ones(1, length(t));
onesM = ones(M, 1);
Mlam = transpose((0:M-1)+lambda);
onesMcotT = onesM*cot(t);
MlamonesT = Mlam*onesT;
twoSinT = onesM*(2*sin(t));
denom = cumprod(twoSinT)./(twoSinT).^(1-lambda);

alpha = onesM*(n*t) + MlamonesT.*(onesM*(t-.5*pi));
if ( ~flag )
    cosAlpha = cos(alpha);
    sinAlpha = sin(alpha);
else
    % Adapted from JACPTS, and LEGPTS:
    %%%%%%%%%%%%%%%% Taylor expansion of cos(alpha0) %%%%%%%%%%%%%%
    k = numel(t):-1:1;
    % HI-LO expansion, to accurately compute (n+.5)*t - (k+.5*lambda-.5)*pi
    ta = double(single(t));
    tb = t - ta;
    hi = n*ta;
    lo = n*tb+lambda*t;
    pia = double(single(pi)); 
    pib = pi - pia;
    dh = (hi-(k-.25)*pia) + lo -.5*(lambda-.5)*pia - (k-.25+.5*(lambda-.5))*pib;

    % Compute cosAlpha(1,:) using Taylor series:
    tmp = 0; sgn = 1; fact = 1; DH = dh; dh2 = dh.*dh;
    for jj = 0:20
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
    for jj = 0:20
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
        cosAlpha(kk,:) = cosAlpha(kk-1,:).*sint + sinAlpha(kk-1,:).*cost;
        sinAlpha(kk,:) = sinAlpha(kk-1,:).*sint - cosAlpha(kk-1,:).*cost;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Sum up all the terms:
vals = C*(c*(cosAlpha./denom));
numer = MlamonesT.*(cosAlpha.*onesMcotT + sinAlpha) + n*sinAlpha;
ders = -C*(c*(numer./denom)); % (dP/dtheta)
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              ASY2 (boundary)                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, w, v, t] = asy2(n, lambda, nbdy, convergence, cte) 
% Algorithm for computing nodes and weights near the boundary.
% N > 2.
% CONVERGENCE ~=0 garantees convergence of Newton-Raphson [4].

k=(1:nbdy).';

% Constant:
Lam = lambda*(1-lambda);

if (convergence ~=0) 
    % This initial guess guarantees convergence [4]:
    t0 = (k-(1-lambda)*0.5)/(n+lambda)*pi;
    cos2 = cos(t0).^2;
    % Sharp initial guess (Förster and Petras, 1993):
    t = t0 + lambda*(1-lambda)/(2*(n+lambda)^2)*(1-(6+lambda*...
        (1-lambda)*(9-2*cos2))./(12*(n+lambda)^2*...
        (1-cos2))).*cot(t0);
    % Choose initial guess that guarantees convergence [4]:
    if  t < t0 
        t = t0;
    end
else
    % Gatteschi's approximation (1979) for initial guesses:
    N = sqrt((n+lambda)^2+lambda*(1-lambda)/3);
    bz = transpose(besselasy(lambda-.5,nbdy)); % Approximates zeros of Bessel
    t = bz/N-Lam/90.*(bz.^3+2*(lambda^2-lambda-.75).*bz)/N^5;
end

% Useful constants:
a=lambda-.5;
b=a;

% Compute higher terms:
[tB1, A2, tB2, A3] = asy2_higherterms(a, b, t, n);

dt = inf; j = 0;
% Newton iteration: 
while ( norm(dt,inf) > sqrt(eps)/1000 && j < 10) 
    [vals, ders] = feval_asy2(n, t, 0); % Evaluate via asymptotic formula.
    dt = vals./ders;                    % Newton update.
    t = t + dt;                         % Next iterate.
    j = j + 1;
end
[vals, ders] = feval_asy2(n, t, 1);     % Evaluate via asymptotic formula.
dt = vals./ders;                        % Newton update
t = t + dt;    
    
% flip:
t = t(nbdy:-1:1); 
ders = ders(nbdy:-1:1);

% Revert to x-space:
x = cos(t);      
w = (1./ders.^2);   

% Compute the constant for weights (adapted from JACPTS):
M = min(20, n-1); 
C1 = 1; 
phi = -a^2/n;
for m = 1:M
    C1 = C1 + phi;
    phi = -(m+a)^2/(m+1)/(n-m)*phi;
    if ( abs(phi/C1) < eps/100 )
        break
    end
end    
C1 = 4^lambda*C1;
w = C1*w;

% Revert to x-space:
v = sqrt(w/cte).*sin(t); 


% The function below is adapted from JACPTS:
function [vals, ders] = feval_asy2(n, t, flag)
% Evaluate asymptotic formula (boundary) - Baratella and Gatteschi (1998). 
    
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
vals2 = Jab + Jbb.*tB0/rho2 + Jab.*A1/rho2^2 + Jbb.*tB1t/rho2^3 + Jab.*A2t/rho2^4;

% Higher terms (not needed for n > 1000).
tB2t = tB2(t); A3t = A3(t);
vals = vals + Jb.*tB2t/rho^5 + Ja.*A3t/rho^6;
vals2 = vals2 + Jbb.*tB2t/rho2^5 + Jab.*A3t/rho2^6;
   
% Constant out the front (Computed accurately!)
ds = .5*(a^2)/n;
s = ds; jj = 1;
while abs(ds/s) > eps/10
    jj = jj+1;
    ds = -(jj-1)/(jj+1)/n*(ds*a);
    s = s + ds;
end
p2 = exp(s)*sqrt((n+a)/n)*(n/rho)^a;
g = [1 1/12 1/288 -139/51840 -571/2488320 163879/209018880 ...
     5246819/75246796800 -534703531/902961561600 ...
     -4483131259/86684309913600 432261921612371/514904800886784000];
f = @(z) sum(g.*[1 cumprod(ones(1,9)./z)]);
C = p2*(f(n+a)/f(n))/sqrt(2);

% Scaling:
valstmp = C*vals;
denom = (sin(t)/2).^lambda;
vals = sqrt(t).*valstmp./denom;

% Relation for derivative:
C2 = C*n/(n+a)*(rho/rho2)^a;
ders = (-n*(2*(n+lambda)-1)*cos(t).*valstmp + 2*(n+a)^2*C2*vals2)/(2*(n+lambda)-1);
ders = ders.*(sqrt(t)./(denom.*sin(t)));
end

end


function [j] = besselasy(v,nbdy) % Approximates zeros of the Besssel 
% function with parameter -1<=V<=5.
% [J] = BESSELASY(V,NBDY) returns the first NBDY approximations.

% Piessens's Chebyshev series approximations (1984). Calculates the 6 first
% zeros to at least 12 decimal figures in region -1<=V<=5:
c1  = [2.883975316228     7.67665211539e-01    -8.6538804759e-02...
    2.0433979038e-02    -6.103761347e-03     2.046841322e-03...
    -7.34476579e-04     2.75336751e-04    -1.06375704e-04...
    4.2003336e-05    -1.6858623e-05     6.85244e-06    -2.8133e-06...
    1.164419e-06    -4.85189e-07     2.03309e-07    -8.5602e-08...
    3.6192e-08    -1.5357e-08     6.537e-09    -2.791e-09     1.194e-09...
    -5.12e-10     2.2e-10    -9.5e-11     4.1e-11    -1.8e-11     8e-12...
    -3e-12     1e-12];

c2 = [8.263194332307    4.209200330779   -0.164644722483...
    0.039764618826   -0.011799527177   0.003893555229  -0.001369989689...
    0.0005030547   -0.00019038177   7.3681222e-05   -2.901083e-05...
    1.1579131e-05  -4.672877e-06  1.903082e-06  -7.8103e-07 3.22648e-07...
    -1.34047e-07 5.5969e-08 -2.3472e-08  9.882e-09 -4.175e-09 1.77e-09...
    -7.52e-10  3.21e-10 -1.37e-10 5.9e-11  -2.5e-11 1.1e-11 -5e-12...
    2e-12];

c3  = [11.493871452173  4.317988625384  -0.130667664397   0.023009510531...
    -0.004987164201   0.001204453026  -0.000310786051  8.383477e-05...
    -2.3343325e-05 6.655551e-06 -1.932603e-06 5.69367e-07  -1.69722e-07...
    5.1084e-08  -1.5501e-08  4.736e-09  -1.456e-09  4.5e-10  -1.4e-10...
    4.3e-11    -1.4e-11  4e-12];

c4 = [14.6890365059305  4.387437455306   -0.109469595763  0.015359574754...
    -0.002655024938   0.000511852711  -0.000105522473  2.2761626e-05...
    -5.071979e-06  1.158094e-06  -2.6948e-07   6.3657e-08   -1.5222e-08...
    3.677e-09    -8.96e-10  2.2e-10 -5.4e-11  1.3e-11  -3e-12 1e-12];

c5 = [17.866882871378   4.435717974422  -0.094492317231 0.011070071951...
    -0.001598668225  0.000257620149  -4.4416219e-05 8.016197e-06...
    -1.495224e-06  2.85903e-07  -5.5734e-08  1.1033e-08 -2.212e-09...
    4.48e-10    -9.2e-11    1.9e-11   -4e-12];

c6 = [21.0347843080875  4.471319438161 -0.083234240394  0.00838807302...
    -0.001042443435  0.000144611721  -2.1469973e-05  3.337753e-06...
    -5.36428e-07   8.8402e-08   -1.4856e-08   2.536e-09  -4.38e-10...
    7.7e-11    -1.4e-11   2e-12];

p = acos((v-2)/3);

cs = cos((0:29)*p);

j(1) = sqrt(v+1)*sum(c1.*cs);
j(2) = sum(c2.*cs);
j(3) = sum(c3.*cs(1:22));
j(4) = sum(c4.*cs(1:20));
j(5) = sum(c5.*cs(1:17));
j(6) = sum(c6.*cs(1:16));


% McMahon's expansion. This expansion gives very accurate approximation 
% for the sth zero (s>=7) in the whole region v>=-1:
s=(7:nbdy);
beta = .25*(2*v+4*s-1)*pi;
mu = 4*v^2;
mu2 = mu*mu;
mu3 = mu2*mu;
mu4 = mu3*mu;
mu5 = mu4*mu;
mu6 = mu5*mu;
mun = mu-1;
A(1,:) = mun./(8*beta);
A(2,:) = mun*(7*mu-31)./(384*beta.^3);
A(3,:) = 4*mun*(83*mu2-982*mu+3779)./(61440*beta.^5);
A(4,:) = 6*mun*(6949*mu3-153855*mu2+1585743*mu-6277237)./(20643840*...
    beta.^7);
A(5,:) = 144*mun*(70197*mu4-2479316*mu3+48010494*mu2-512062548*mu+...
    2092163573)./(11890851840*beta.^9);
A(6,:) = 720*mun*(5592657*mu5-287149133*mu4+8903961290*mu3-179289628602*...
    mu2+1982611456181*mu-8249725736393)./(10463949619200*beta.^11);
A(7,:) = 576*mun*(4148944183*mu6-291245357370*mu5+13172003634537*mu4...
    -426353946885548*mu3+8929489333108377*mu2-100847472093088506*mu...
    +423748443625564327)./(13059009124761600*beta.^13);

j(s) = beta - sum(A,1);
end


% This code is from JACPTS:
function Ja = besselTaylor(t, z, a)
% BESSELTAYLOR    Accurate evaluation of Bessel function J_A for asy2. (See [2].)
% BESSELTAYLOR(T, Z, A) evaluates J_A(Z+T) by a Taylor series expansion about Z. 

npts = numel(t);
kmax = min(ceil(abs(log(eps)/log(norm(t, inf)))), 30);
H = bsxfun(@power, t, 0:kmax).';
% Compute coeffs in Taylor expansions about z (See NIST 10.6.7)
[nu, JK] = meshgrid(-kmax:kmax, z);
Bjk = besselj(a + nu, JK, 0);
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


% This code is from JACPTS:
function [tB1, A2, tB2, A3, tB3, A4] = asy2_higherterms(alph, bet, theta, n)
% ASY2_HIGHERTERMS   Higher-order terms for boundary asymptotic series.
% Compute the higher order terms in asy2 boundary formula. See [2]. 

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

if ( nargout < 3 )
    % Make function for output
    tB1 = @(theta) bary(theta,tB1,t,v);
    A2 = @(theta) bary(theta,A2,t,v);
    return
end

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

if ( nargout < 6 )
    % Make function for output
    tB1 = @(theta) bary(theta, tB1, t, v);
    A2 = @(theta) bary(theta, A2, t, v);
    tB2 = @(theta) bary(theta, tB2, t, v);
    A3 = @(theta) bary(theta, A3, t, v);
    return
end

% A2p:
A3p = D*A3;
A3p = A3p - A3p(1);
A3p_t = A3p./t;
% Extrapolate point at t = 0:
w = pi/2-t(2:end);
w(2:2:end) = -w(2:2:end);
A3p_t(1) = sum(w.*A3p_t(2:end))/sum(w);

% B2:
tB3 = -.5*A3p - (.5+alph)*(C*A3p_t) + .5*C*(f.*A3);
B3 = tB3./t;
% Extrapolate point at t = 0
B3(1) = sum(w.*B3(2:end))/sum(w);

% A3:
K = C*(f.*tB3);
A4 = .5*(D*tB3) - (.5+alph)*B3 - .5*K;
A4 = A4 - A4(1);

% Make function for output:
tB1 = @(theta) bary(theta, tB1, t, v);
A2 = @(theta) bary(theta, A2, t, v);
tB2 = @(theta) bary(theta, tB2, t, v);
A3 = @(theta) bary(theta, A3, t, v);
tB3 = @(theta) bary(theta, tB3, t, v);
A4 = @(theta) bary(theta, A4, t, v);

end
