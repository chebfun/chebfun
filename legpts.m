function [x, w, v, t] = legpts(n, int, meth)
%LEGPTS    Legendre points and Gauss-Legendre Quadrature Weights.
%   LEGPTS(N) returns N Legendre points X in (-1,1).
%
%   [X, W] = LEGPTS(N) returns also a row vector W of weights for Gauss-Legendre
%   quadrature.
%
%   [X, W] = LEGPTS(N, INTERVAL) scales the nodes and weights for the finite
%   interval INTERVAL.
%
%   [X, W, V] = LEGPTS(N) or [X, W, V] = LEGPTS(N, D) returns additionally a
%   column vector V of weights in the barycentric formula corresponding to the
%   points X. The weights are scaled so that max(abs(V)) = 1.
%
%   LEGPTS(N, INTERVAL, METHOD) or LEGPTS(N, METHOD) allows the user to select
%   which method to use.
%    METHOD = 'REC' uses the recurrence relation for the Legendre polynomials
%     and their derivatives to perform Newton iteration on the WKB approximation
%     to the roots. Default for N < 100.
%    METHOD = 'ASY' uses the Hale-Townsend fast algorithm based upon asymptotic
%     formulae, which is fast and accurate. Default for N >= 100.
%    METHOD = 'GW' will use the traditional Golub-Welsch eigenvalue method,
%       which is maintained mostly for historical reasons.
%
%   [X, W, V, T] = LEGPTS(N) returns also the arccos of the nodes, T = acos(X).
%   In some situations (in particular with 'ASY') these can be computed to a
%   much better relative precision than X.
%
% See also CHEBPTS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  'GW' by Nick Trefethen, March 2009 - algorithm adapted from [1].
%  'REC' by Nick Hale, July 2011
%  'ASY' by Nick Hale & Alex Townsend, May 2012 - see [2].
%
%  References:
%   [1] G. H. Golub and J. A. Welsch, "Calculation of Gauss quadrature rules",
%       Math. Comp. 23:221-230, 1969,
%   [2] N. Hale and A. Townsend, "Fast computation of Gauss-Jacobi quadrature
%       nodes and weights",In preparation, 2012.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Defaults
interval = [-1, 1];
method = 'default';
method_set = nargin == 3;

% Deal with trivial cases:
if ( n < 0 )
    error('CHEBFUN:legpts:nNegative', ...
        'First input should be a positive number.');
elseif ( n == 0 )   % Return empty vectors if n == 0:
    x = [];
    w = [];
    v = [];
    t = [];
    return
elseif ( n == 1 )
    x = 0;
    w = 2;
    v = 1;
    t = 1;
    return
elseif ( n == 2 )
    x = [-1 ; 1]/sqrt(3);
    w = [1 1];
    v = [1 ; -1];
    t = acos(x);
    return
end

% Check the inputs:
if ( nargin > 1 )
    if ( nargin == 3 )
        % Calling sequence = LEGPTS(N, INTERVAL, METHOD)
        interval = int;
        method = meth;
    elseif ( nargin == 2 )
        if ( ischar(int) )
            % Calling sequence = LEGPTS(N, METHOD)
            method = int;
            method_set = true;
        else
            % Calling sequence = LEGPTS(N, INTERVAL)
            interval = int;
        end
    end
    validStrings = {'default', 'GW', 'ASY', 'REC'};
    if ( ~any(strcmpi(method, validStrings)) )
        if ( strcmpi(method, 'GLR') )
            error('CHEBFUN:legpts:glr', ...
                'The GLR algorithm is no longer supported.');
        end
        error('CHEBFUN:legpts:inputs', ['Unrecognised input string: ', method]);
    end
    if ( numel(interval) > 2 )
        warning('CHEBFUN:legpts:domain',...
            'Piecewise intervals are not supported and will be ignored.');
        interval = interval([1, end]);
    end
end
if ( any(isinf(interval)) )
    error('CHEBFUN:legpts:interval', 'Unbounded intervals are not supported.');
end

if ( n <= 20 )
    % Force REC for n <= 20:
    method = 'rec'; % 
    method_set = 1; 
end

% Choose the method:
t = [];
if ( (n < 100 && ~method_set) || strcmpi(method, 'rec') )
    [x, w, v] = rec(n);        % REC (Standard recurrence relation)
elseif ( strcmpi(method, 'GW') )
    [x, w, v] = gw(n);         % GW see [1]
else
    [x, w, v, t] = asy(n);     % HT see [2]
end

% Normalise the barycentric weights:
v = abs(v);
v = v./max(v);
v(2:2:end) = -v(2:2:end);

% Compute a T is one is asked for:
if ( nargout == 4 && isempty(t) )
    t = acos(x);
end

% Rescale to arbitrary finite interval:
if ( ~all(interval == [-1 1]) )
    dab = diff(interval);
    x = (x+1)/2*dab + interval(1);
    w = dab*w/2;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------------ Routines for GW algorithm ---------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, w, v] = gw(n)
beta = .5./sqrt(1-(2*(1:n-1)).^(-2)); % 3-term recurrence coeffs
T = diag(beta, 1) + diag(beta, -1);   % Jacobi matrix
[V,D] = eig(T);                       % Eigenvalue decomposition
x = diag(D);                          % Legendre points
[x, i] = sort(x);                     % Sort
w = 2*V(1,i).^2;                      % Quadrature weights
v = sqrt(1-x.^2).*abs(V(1,i))';       % Barycentric weights

% Enforce symmetry:
ii = 1:floor(n/2);
x = x(ii);
w = w(ii);
vmid = v(floor(n/2) + 1);
v = v(ii);
if ( mod(n, 2) )
    % Odd number.
    x = [x ; 0 ; -x(end:-1:1)];
    w = [w,  2 - sum(2*w), w(end:-1:1)];
    v = [v ; vmid ; v(end:-1:1)];
else
    % Evem number.
    x = [x ; -x(end:-1:1)];
    w = [w, w(end:-1:1)];
    v = [v ; v(end:-1:1)];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ----------------------- Routines for REC algorithm ---------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, w, v] = rec(n)

% Asymptotic formula (Tricomi) - only for positive x.
if ( mod(n,2) )
    s = 1;
else
    s = 0;
end
k = ((n+s)/2:-1:1).';
theta = pi*(4*k-1)/(4*n+2);
x = ( 1 - (n-1)/(8*n^3) - 1/(384*n^4)*(39-28./sin(theta).^2) ).*cos(theta);

% Initialise:
Pm2 = 1;
Pm1 = x;
PPm2 = 0;
PPm1 = 1;
dx = inf;
counter = 0;

% Loop until convergence:
while ( norm(dx, inf) > eps && counter < 10 )
    counter = counter + 1;
    for k = 1:n-1,
        P = ((2*k+1)*Pm1.*x-k*Pm2)/(k+1);
        Pm2 = Pm1;
        Pm1 = P;
        PP = ((2*k+1)*(Pm2+x.*PPm1)-k*PPm2)/(k+1);
        PPm2 = PPm1;
        PPm1 = PP;
    end
    % Newton step:
    dx = -P./PP;
    % Newton update:
    x = x + dx;
    % Reinitialise:
    Pm2 = 1;
    Pm1 = x;
    PPm2 = 0;
    PPm1 = 1;
end

% Once more for derivatives:
for k = 1:n-1,
    P = ( (2*k+1)*Pm1.*x - k*Pm2 ) / (k+1);
    Pm2 = Pm1;
    Pm1 = P;
    PP = ( (2*k+1)*(Pm2+x.*PPm1) - k*PPm2 ) / (k+1);
    PPm2 = PPm1;
    PPm1 = PP;
end

% [TOD0]: This relation might prove useful?
%     PP = -n*(x.*P-Pm2)./(1-x.^2);

% Reflect for negative values:
x = [-x(end:-1:1+s) ; x];
ders = [PP(end:-1:1+s) ; PP];

% Quadrature weights:
w = 2./((1-x.^2).*ders.^2)';

% Barycentric weights
v = 1./ders;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------------- Routines for ASY algorithm ------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, w, v, t] = asy(n)

% Determine switch between interior and boundary regions:
nbdy = min(10, floor(n/2)); % Typically, the 10 nodes nearest the boundary.

% Interior algorithm:
[x, w, v, t] = asy1(n, nbdy);

% Boundary algorithm:
[xbdy, wbdy, vbdy, tbdy] = asy2(n, nbdy);

% Combine:
bdyidx1 = n-(nbdy-1):n;
bdyidx2 = nbdy:-1:1;
x(bdyidx1) = xbdy;
w(bdyidx1) = wbdy;
v(bdyidx1) = vbdy;
t(bdyidx1) = tbdy;

% Reflect using symmetry:
x(bdyidx2) = -xbdy;
w(bdyidx2) = wbdy;
v(bdyidx2) = vbdy;
t(bdyidx2) = -tbdy;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             ASY1 (Interior)                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, w, v, t] = asy1(n, nbdy)
%ASY1 Interior asymptotics method.

%%%%%%%%%%%%%% Approximate roots via asymptotic formula. %%%%%%%%%%
% Tricomi's formula:
modn2 = mod(n, 2);
k = (n-2+modn2)/2+1:-1:1;
theta = pi*(4*k-1)/(4*n+2);
x = ( 1 - (n-1)/(8*n^3) - 1/(384*n^4)*(39-28./sin(theta).^2) ).*cos(theta);
t = acos(x);   % Approximate roots in theta-space.
if ( n < 666 ) % <-- Determined experimentally.
    % If n is small, use Olver's (1974) apprioximation in the interior.
    idx = (x > .5);
    npts = sum(idx);
    % Roots of the Bessel function J_0 (Precomputed in Mathematica)
    jk = [2.404825557695773     5.520078110286310    8.653727912911012 ...
        11.791534439014281    14.930917708487785   18.071063967910922 ...
        21.211636629879258    24.352471530749302   27.493479132040254 ...
        30.634606468431975    33.775820213573568].';
    if ( npts > 11 )
        % Estimate the larger Bessel roots: (See Branders et al., JCP 1981).
        p = ( (length(jk)+1:npts).' -.25 )*pi;
        pp = p.*p;
        num = 0.0682894897349453 + pp.*(0.131420807470708 + ...
            pp.*(0.0245988241803681 + pp.*0.000813005721543268));
        den = p.*(1.0 + pp.*(1.16837242570470 + pp.*(0.200991122197811 + ...
            pp.*(0.00650404577261471))));
        jk = [jk ; p + num./den];
    end
    phik = jk(1:npts)/(n+.5);
    tnew = phik + (phik.*cot(phik) - 1)./(8*phik*(n+.5)^2);
    t(idx) = tnew(end:-1:1); % Approximate roots in theta-space.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Locate the boundary node:
mint = t(end - nbdy + 1);
idx = max(find(t < mint, 1) - 1, 1);

% Initialise:
dt = inf;
% Newton iteration: (Always converges)
while ( norm(dt, inf) > sqrt(eps)/1000 )       % <-- Enough, as again below
    [vals, ders] = feval_asy1(n, t , mint, 0); % Evaluate via asy formulae
    dt = vals./ders;                           % Newton step
    t = t - dt;                                % Newton update
    dt = dt(1:idx-1);                          % Ignore boundary terms
end
[vals, ders] = feval_asy1(n, t, mint, 1);  % Once more for good ders.
t = t - vals./ders;                        % Newton update

% Compute x, w, and v:
x = cos(t);
w = 2./ders.^2;
v = sin(t)./ders;

% Flip using symmetry for negative nodes:
if ( modn2 )
    x = [-x(end:-1:2), x].';
    w = [w(end:-1:2), w];
    v = -[v(end:-1:2), v].';
    t = [-t(end:-1:2), t].';
else
    x = [-x(end:-1:1), x].';
    w = [w(end:-1:1), w];
    v = [-v(end:-1:1), v].';
    t = [-t(end:-1:1), t].';
end

end

function [vals, ders] = feval_asy1(n, t , mint, flag) %#ok<INUSD>
%FEVAL_ASY1   Evaluate 1st asymptotic formula (interior)

% Max number of expansion terms:
M = 20; % <-- Provides 16 digit accuracy for n >= 100.

% Coefficients in expansion:
c = cumprod( (1:2:2*M-1)./(2:2:2*M) );
d = cumprod( (1:2:2*M-1)./(2*n+3:2:2*(n+M)+1) );
c = [1, c.*d];

% How many terms required in the expansion?
R = (8/pi)*c./(2*sin(mint)).^(.5:M+1)/10;
R = R(abs(R) > eps);
M = length(R);
c = c(1:M);

%%%%%% Constant out the front: ( C = sqrt(4/pi)*gamma(n+1)/gamma(n+3/2) ) %%%%
ds = -1/8/n;
s = ds;
j = 1;
while ( abs(ds/s) > eps/100 ) % Taylor series in expansion 
    j = j + 1;
    ds = -.5*(j-1)/(j+1)/n*ds;
    s = s + ds;
end
p2 = exp(s)*sqrt(4/(n+.5)/pi);
% Stirling's series:
g = [1, 1/12, 1/288, -139/51840, -571/2488320, 163879/209018880, ...
    5246819/75246796800, -534703531/902961561600, ...
    -4483131259/86684309913600, 432261921612371/514904800886784000];
f = @(z) sum(g.*[1, cumprod(ones(1, 9)./z)]);
C = p2*(f(n)/f(n+.5));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Some often used vectors/matrices:
onesT = ones(1, length(t));
onesM = ones(M, 1);
M05 = transpose((0:M-1) +.5);
onesMcotT = onesM*cot(t);
M05onesT = M05*onesT;
twoSinT = onesM*(2*sin(t));
denom = cumprod(twoSinT)./sqrt(twoSinT);

alpha = onesM*(n*t) + M05onesT.*(onesM*(t-.5*pi));
if ( ~flag )
    cosAlpha = cos(alpha);
    sinAlpha = sin(alpha);
else
    %%%%%%%%%%%%%%%% Taylor expansion of cos(alpha0) %%%%%%%%%%%%%%
    k = numel(t):-1:1;
    rho = n + 0.5;
    % HI-LO expansion, to accurately compute rho*t - (k-.25)*pi
    ta = double(single(t));
    tb = t - ta;
    hi = rho*ta;
    lo = rho*tb;
    pia = double(single(pi));
    pib = -8.742278000372485e-08; %pib = pi - pia;
    dh = (hi-(k-.25)*pia) + lo - (k-.25)*pib;

    % Compute cosAlpha(1,:) using Taylor series:
    tmp = 0; sgn = 1; fact = 1; DH = dh; dh2 = dh.*dh;
    for j = 0:20
        dc = sgn*DH/fact;
        tmp = tmp + dc;
        sgn = -sgn;
        fact = fact*(2*j+3)*(2*j+2);
        DH = DH.*dh2;
        if ( norm(dc,inf) ) < eps/2000, break, end
    end
    tmp(2:2:end) = -tmp(2:2:end);
    tmp = sign(cos((n+.5)*t(2)-.25*pi)*tmp(2))*tmp;
    cosAlpha(1,:) = tmp;

    % Compute sinAlpha(1,:) using Taylor series:
    tmp = 0; sgn = 1; fact = 1; DH = 1; dh2 = dh.*dh;
    for j = 0:20
        dc = sgn*DH/fact;
        tmp = tmp + dc;
        sgn = -sgn;
        fact = fact*(2*j+2)*(2*j+1);
        DH = DH.*dh2;
        if (norm(dc, inf)) < eps/2000, break, end
    end
    tmp(2:2:end) = -tmp(2:2:end);
    tmp = sign(sin((n+.5)*t(2)-.25*pi)*tmp(2))*tmp;
    sinAlpha(1,:) = tmp;

    % Compute cosAlpha(k,:) and sinAlpha(k,:) for k = 2,...,M:
    sint = sin(t);
    cost = cos(t);
    for k = 2:M
        cosAlpha(k,:) = cosAlpha(k-1,:).*sint+sinAlpha(k-1,:).*cost;
        sinAlpha(k,:) = sinAlpha(k-1,:).*sint-cosAlpha(k-1,:).*cost;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% Sum up all the terms:
vals = C*(c*(cosAlpha./denom));
numer = M05onesT.*(cosAlpha.*onesMcotT + sinAlpha) + n*sinAlpha;
ders = -C*(c*(numer./denom)); % (dP/dtheta)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             ASY2 (Boundary)                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, w, v, t, ders] = asy2(n, npts)
%ASY2 Boundary asymptotics method.

if ( npts > ceil((n+1)/2) )
    error('CHEBFUN:legpts:N', 'NPTS must be <= N/2');
end

%%%%%%%% Approximation for Legendre roots: (See Olver 1974) %%%%%%%%
% Roots pf the Bessel function J_0: (Precomputed in Mathematica)
jk = [2.404825557695773,     5.520078110286310,    8.653727912911012, ...
    11.791534439014281,    14.930917708487785,   18.071063967910922, ...
    21.211636629879258,    24.352471530749302,   27.493479132040254, ...
    30.634606468431975,    33.775820213573568].';
phik = jk(1:npts)/(n+.5);
t = phik + (phik.*cot(phik)-1)./(8*phik*(n+.5)^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[tB1, A2, tB2, A3] = asy2_higherterms(0, 0, t, n);
dt = inf;
% Newton iteration: (Always converges)
while ( norm(dt,inf) > sqrt(eps)/1000 )  % <-- Enough as once more below
    [vals, ders] = feval_asy2(n, t, 0);  % Evaluate via asy formula
    dt = vals./ders;                     % Newton update
    t = t + dt;                          % Next iterate
end
% Once more for good ders.
[vals, ders] = feval_asy2(n, t, 1);

% Flip:
t = t(npts:-1:1);
ders = ders(npts:-1:1);

% Revert to x-space:
x = cos(t);
w = (2./ders.^2).';
v = sin(t)./ders;

    function [vals, ders] = feval_asy2(n, t, flag)
        %FEVAL_ASY2   Evaluate 2nd asymptotic formula (boundary)
        % If FLAG == 1, then a more accurate formula is used to evaluate J0(r*t)
        
        % Useful constants:
        rho = n + .5;
        rho2 = n - .5;
        
        % Evaluate the Bessel functions:
        Ja = besselj(0, rho*t, 0);
        Jb = besselj(1, rho*t, 0);
        Jbb = besselj(1, rho2*t, 0);
        if ( ~flag )
            Jab = besselj(0,rho2*t,0);
        else
            % In the final step, perform accurate evaluation:
            Jab = besselTaylor(-t, rho*t);
        end
        
        % Evaluate functions for recurrsive definition of coefficients:
        gt = .5*(cot(t) - 1./t);
        gtdt = .5*(-csc(t).^2 + 1./t.^2);
        tB0 = .25*gt;
        A1 = gtdt/8 - 1/8*gt./t - gt.^2/32;
        % Higher terms
        tB1t = tB1(t);
        A2t = A2(t);
        
        % Values:
        vals = Ja + Jb.*tB0/rho + Ja.*A1/rho^2 + Jb.*tB1t/rho^3 + Ja.*A2t/rho^4;
        % derivatives:
        vals2 = Jab + Jbb.*tB0/rho2 + Jab.*A1/rho2^2 + ...
            Jbb.*tB1t/rho2^3 + Jab.*A2t/rho2^4;
        
        % Higher terms: (not needed for n > 1000)
        tB2t = tB2(t);
        A3t = A3(t);
        vals = vals + Jb.*tB2t/rho^5 + Ja.*A3t/rho^6;
        vals2 = vals2 + Jbb.*tB2t/rho2^5 + Jab.*A3t/rho2^6;
        
        % Relation for derivative:
        ders = n*(-cos(t).*vals + vals2)./sin(t);
        
        % Common factors:
        denom = sqrt(t./sin(t));
        ders = ders.*denom;
        vals = vals.*denom;
        
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [TODO]: The following are duplicated in LEGPTS() and JACPTS().
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Ja = besselTaylor(t, z, a)
%BESSELTAYLOR    Accurate evaluation of Bessel function J0 for asy2. (See [2].)
% BESSELTAYLOR(T, Z) evaluates J0(Z+T) by a Taylor series expansion about Z. 

npts = numel(t);
if ( nargin < 3 ), a = 0; end
kmax = min(ceil(abs(log(eps)/log(norm(t, inf)))), 30);
H = bsxfun(@power, t, 0:kmax).';
% Compute coeffs in Taylor expansions about z: (See NIST 10.6.7)
[nu, JK] = meshgrid(-kmax:kmax, z);
Bjk = besselj(a + nu, JK, 0);
nck = abs(pascal(floor(1.25*kmax), 1)); nck(1,:) = []; % nchoosek
AA = [Bjk(:,kmax+1),  zeros(npts,kmax)];
kFactorial = 1;
for k = 1:kmax
    sgn = 1;
    for l = 0:k
        AA(:,k+1) = AA(:,k+1) + sgn*nck(k,l+1)*Bjk(:,kmax+2*l-k+1);
        sgn = -sgn;
    end
    kFactorial = k*kFactorial;
    AA(:,k+1) = AA(:,k+1)/2^k/kFactorial;
end
% Evaluate Taylor series:
Ja = zeros(npts,1);
for k = 1:npts
    Ja(k,1) = AA(k,:)*H(:,k);
end

end

function [tB1, A2, tB2, A3, tB3, A4] = asy2_higherterms(alph, bet, theta, n)
%ASY2_HIGHERTERMS   Higher-order terms for boundary asymptotic series.
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
Nm1 = N-1;

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
C = colloc2.cumsummat(N)*(.5*c);
D = colloc2.diffmat(N)*(2/c);
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
    tB1 = @(theta) chebtech.bary(theta,tB1,t,v);
    A2 = @(theta) chebtech.bary(theta,A2,t,v);
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
    tB1 = @(theta) chebtech.bary(theta, tB1, t, v);
    A2 = @(theta) chebtech.bary(theta, A2, t, v);
    tB2 = @(theta) chebtech.bary(theta, tB2, t, v);
    A3 = @(theta) chebtech.bary(theta, A3, t, v);
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
tB1 = @(theta) chebtech.bary(theta, tB1, t, v);
A2 = @(theta) chebtech.bary(theta, A2, t, v);
tB2 = @(theta) chebtech.bary(theta, tB2, t, v);
A3 = @(theta) chebtech.bary(theta, A3, t, v);
tB3 = @(theta) chebtech.bary(theta, tB3, t, v);
A4 = @(theta) chebtech.bary(theta, A4, t, v);

end
