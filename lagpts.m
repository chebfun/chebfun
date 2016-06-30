function [x, w, v] = lagpts(n, int, meth)
%LAGPTS  Laguerre points and Gauss-Laguerre Quadrature Weights.
%   LAGPTS(N) returns N Laguerre points X in (0,inf).
%
%   [X, W] = LAGPTS(N) returns also a row vector W of weights for Gauss-Laguerre
%   quadrature. [X, W, V] = LAGPTS(N) returns in addition a column vector V
%   of the barycentric weights corresponding to X.
%
%   LAGPTS(N, D) scales the nodes and weights for the semi-infinite domain D.
%   D can be either a domain object or a vector with two components.
%
%   [X, W] = LAGPTS(N, METHOD) allows the user to select which method to use.
%   METHOD = 'GW' will use the traditional Golub-Welsch eigenvalue method,
%   which is best for when N is small. METHOD = 'FAST' will use the
%   Glaser-Liu-Rokhlin fast algorithm, which is much faster for large N.
%   METHOD = 'RH' will use asymptotics of Laguerre polynomials, and METHOD =
%   'RHW' is O(sqrt(n)) as it stops when the weights are below realmin.
%   By default LAGPTS uses 'GW' when N < 128, 'FAST' when N is in [128, 4200]
%   and else 'RH'.
%
% References:
%   [1] G. H. Golub and J. A. Welsch, "Calculation of Gauss quadrature rules",
%       Math. Comp. 23:221-230, 1969,
%   [2] A. Glaser, X. Liu and V. Rokhlin, "A fast algorithm for the calculation
%       of the roots of special functions", SIAM Journal on Scientific
%       Computing", 29(4):1420-1438, 2007.
%   [3] P. Opsomer, (in preparation).
%   [4] M. Vanlessen, "Strong asymptotics of Laguerre-Type orthogonal
%       polynomials and applications in Random Matrix Theory", Constr. Approx.,
%       25:125-175, 2007.
%
% See also CHEBPTS, LEGPTS, HERMPTS, and JACPTS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
%
% 'GW' by Nick Trefethen, March 2009 - algorithm adapted from [1].
% 'FAST' by Nick Hale, March 2010 - algorithm adapted from [2].
% 'RH' by Peter Opsomer, June 2016 - algorithm adapted from [3], based on [4].

% Defaults:
method = 'default';
interval = [0, inf];

if ( n < 0 )
    error('CHEBFUN:lagpts:n', ...
        'First input should be a positive integer.');
end

% Return empty vector if n = 0.
if ( n == 0 )
    [x, w, v] = deal([]);
    return
end

% Check the inputs.
if ( nargin > 1 )
    if ( nargin == 3 )
        interval = int;
        method = meth;
    elseif ( nargin == 2 )
        if ( ischar(int) )
            method = int;
        else
            interval = int;
        end
    end
    if ( ~any(strcmpi(method, {'default', 'GW', 'fast', 'RH', 'RHW'})) )
        error('CHEBFUN:lagpts:inputs', 'Unrecognised input string %s.', method);
    end
    if ( numel(interval) > 2 )
        warning('CHEBFUN:lagpts:domain',...
            'Piecewise intervals not supported and will be ignored.');
        interval = interval([1, end]);
    end
end

if ( sum(isinf(interval)) ~= 1 )
    error('CHEBFUN:lagpts:inf', 'LAGPTS only supports semi-infinite domains.');
end

% Decide to use GW, FAST or RH.
if ( strcmpi(method,'GW') || ( ( n < 128 ) && strcmpi(method,'default') ) )
    % GW, see [1]
    alpha = 2*(1:n)-1;  beta = 1:n-1;     % 3-term recurrence coeffs
    T = diag(beta,1) + diag(alpha) + diag(beta,-1);  % Jacobi matrix
    [V, D] = eig(T);                      % eigenvalue decomposition
    [x, indx] = sort(diag(D));            % Laguerre points
    w = V(1,indx).^2;                     % Quadrature weights
    v = sqrt(x).*abs(V(1,indx)).';        % Barycentric weights
    v = v./max(v);
    v(2:2:n) = -v(2:2:n);
    
elseif ( strcmpi(method,'fast') || ( ( n < 4200) && strcmpi(method,'default') ) )
    % Fast, see [2]
    [x, ders] = alg0_Lag(n);              % Nodes and L_n'(x)
    w = exp(-x)./(x.*ders.^2); w = w';    % Quadrature weights
    v = exp(-x/2)./ders;                  % Barycentric weights
    v = -v./max(abs(v));
    
else
    % RH, see [3] and [4]
    [x, w] = alg_rh(n, strcmpi(method, 'RHW') );  % Nodes and quadrature weights
    v = sqrt(w'.*x);                              % Barycentric weights
    v = -v./max(abs(v));
    
end
w = (1/sum(w))*w;                                 % Normalise so that sum(w) = 1

% Nonstandard interval
if ( ~all(interval == [0, inf]) )
    a = interval(1);
    b = interval(2);
    if ( isinf(b) )
        x = x + a;
        w = w*exp(-a);
    else
        x = -x + b;
        w = w*exp(b);
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%% Routines for FAST algorithm %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, ders] = alg0_Lag(n)
ders = zeros(n, 1);
xs = 1/(2*n+1);
n1 = 20;
n1 = min(n1, n);
x = zeros(n, 1);
for k = 1:n1
    [xs, ders(k)] = alg3_Lag(n, xs);
    x(k) = xs;
    xs = 1.1*xs;
end
[x, ders] = alg1_Lag(x, ders, n, n1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [roots, ders] = alg1_Lag(roots, ders, n, n1)

% number of terms in Taylor expansion
m = 30;

% initialise
hh1 = ones(m+1, 1);
zz = zeros(m, 1);
u = zeros(1, m+1);
up = zeros(1, m+1);

x = roots(n1);
for j = n1:(n - 1)
    
    % initial approx
    h = rk2_Lag(pi/2, -pi/2, x, n) - x;
    
    % scaling:
    M = 1/h;
    M2 = M^2;
    M3 = M^3;
    M4 = M^4;
    
    % recurrence relation for Laguerre polynomials
    r = x*(n + .5 - .25*x);
    p = x^2;
    u(1:2) = [0 ; ders(j)/M];
    u(3) = -.5*u(2)/(M*x) - (n + .5 - .25*x)*u(1)/(x*M^2);
    u(4) = -u(3)/(M*x) + ( -(1+r)*u(2)/6/M^2 - (n+.5-.5*x)*u(1)/M^3 ) / p;
    up(1:3) = [u(2) ; 2*u(3)*M ; 3*u(4)*M];
    
    for k = 2:(m - 2)
        u(k+3) = ( -x*(2*k+1)*(k+1)*u(k+2)/M - (k*k+r)*u(k+1)/M2 - ...
            (n+.5-.5*x)*u(k)/M3 + .25*u(k-1)/M4 ) / (p*(k+2)*(k+1));
        up(k+2) = (k+2)*u(k+3)*M;
    end
    up(m+1) = 0;
    
    % Flip for more accuracy in inner product calculation.
    u = u(m+1:-1:1);
    up = up(m+1:-1:1);
    
    % Newton iteration
    hh = hh1;
    hh(end) = M;
    step = inf;
    l = 0;
    if ( M == 1 )
        Mhzz = (M*h) + zz;
        hh = [M ; cumprod(Mhzz)];
        hh = hh(end:-1:1);
    end
    while ( (abs(step) > eps) && (l < 10) )
        l = l + 1;
        step = (u*hh)/(up*hh);
        h = h - step;
        Mhzz = (M*h) + zz;
        % Powers of h (This is the fastest way!)
        hh = [M ; cumprod(Mhzz)];
        % Flip for more accuracy in inner product
        hh = hh(end:-1:1);
    end
    
    % Update
    x = x + h;
    roots(j+1) = x;
    ders(j+1) = up*hh;
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x1, d1] = alg3_Lag(n, xs)
[u, up] = eval_Lag(xs, n);
theta = atan(sqrt(xs/(n + .5 - .25*xs))*up/u);
x1 = rk2_Lag(theta, -pi/2, xs, n);

% Newton iteration
step = inf;
l = 0;
while ( (abs(step) > eps || abs(u) > eps) && (l < 200) )
    l = l + 1;
    [u, up] = eval_Lag(x1, n);
    step = u/up;
    x1 = x1 - step;
end

[~, d1] = eval_Lag(x1, n);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Evaluate Laguerre polynomial via recurrence
function [L, Lp] = eval_Lag(x, n)
Lm2 = 0;
Lm1 = exp(-x/2);
Lpm2 = 0;
Lpm1 = 0;
for k = 0:n-1
    L = ( (2*k+1-x).*Lm1 - k*Lm2 ) / (k + 1);
    Lp = ( (2*k+1-x).*Lpm1 - Lm1 - k*Lpm2 ) / (k + 1);
    Lm2 = Lm1;
    Lm1 = L;
    Lpm2 = Lpm1;
    Lpm1 = Lp;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Runge-Kutta for Laguerre Equation
function x = rk2_Lag(t, tn, x, n)
m = 10;
h = (tn - t)/m;
for j = 1:m
    f1 = (n + .5 - .25*x);
    k1 = -h/( sqrt(f1/x) + .25*(1/x-.25/f1)*sin(2*t) );
    t = t + h;
    x = x + k1;
    f1 = (n + .5 - .25*x);
    k2 = -h/( sqrt(f1/x) + .25*(1/x-.25/f1)*sin(2*t) );
    x = x + .5*(k2 - k1);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%% Routines for RH algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, w] = alg_rh(n, compRepr)

if compRepr
    % Get a heuristic for the indices where the weights are about above realmin.
    mn = min(n,ceil(17*sqrt(n)));
else
    mn = n;
end
% Initial guesses
x = [ besselroots(0, 3).^2/(4*n + 2); zeros(mn-3, 1) ];
w = zeros(1, mn);

factor0 = 2.757254379232566e-04*n^(-6) + 1.511212766999818e-03*n^(-5) - ...
    8.937757201646138e-04*n^(-4) - 4.783950617283954e-03*n^(-3) +...
    1.388888888888890e-02*n^(-2) + 1.666666666666667e-01*n^(-1) + 1;
factor1 = 1.786938204923081e-03*(n-1)^(-6) + 6.174370468351152e-04*(n-1)^(-5)...
    - 5.677726337448679e-03*(n-1)^(-4) + 9.104938271604964e-03*(n-1)^(-3) + ...
    1.805555555555556e-01*(n-1)^(-2) + 1.166666666666667*(n-1)^(-1) + 1;
% We factored out some constants from the ratio or product of the asymptotics.
factorx = sqrt(factor1/factor0)/(2 - 2/n);
factorw = -(1 - 1/(n + 1) )^(n + 1)*(1 - 1/n)*exp(1 + 2*log(2) )*4*pi*...
    sqrt(factor0*factor1);

% This is a heuristic for the number of terms in the expansions that follow.
T = ceil(25/log(n) );
% Start with the expansion in terms of Bessel functions. We could also define
% poly as a function which switches between pl, pb and pr, or triple the loop.
poly = @pl;
for k = 1:mn
    if ( k > 3 ) % Use quadratic extrapolation for the initial guesses.
        x(k) = 3*x(k-1) - 3*x(k-2) + x(k-3);
    end
    if x(k) > 3.7*n
        poly = @pr;
    elseif x(k) > sqrt(n)
        % The fixed delta in the RHP would mean this bound is proportional to
        % n, but x(1:k) are O(1/n) so choose bound in between to
        % make more use of the (cheaper) expansion in the bulk.
        poly = @pb;
    end
    step = x(k);
    l = 0; % Newton-Raphson iteration number
    ov = inf; % Previous/old value
    ox = x(k); % Old x
    % [FIXME] Accuracy of the expansions up to eps would lower this bound.
    while ( ( abs(step) > eps*400*x(k) ) && ( l < 20) )
        l = l + 1;
        pe = poly(n, x(k), 0, T);
        % poly' = (p*exp(-Q/2) )' = exp(-Q/2)*(p' -p/2) with orthonormal p
        step = pe/(poly(n-1, x(k), 1, T)*factorx - pe/2);
        if (abs(pe) >= abs(ov)*(1-5e5*eps) )
            % The function values do not decrease enough any more due to
            % roundoff errors, so set to the previous value and quit.
            x(k) = ox;
            break
        end
        ox = x(k);
        x(k) = x(k) -step;
        ov = pe;
    end
    w(k) = factorw/poly(n - 1, x(k), 1, T)/poly(n+1, x(k), 0, T)/exp( x(k) );
    if ( w(k) == 0 ) && ( k > 1 ) && ( w(k-1) > 0 ) % We could stop now.
        if compRepr
            w = w(1:k-1);
            x = x(1:k-1);
            return;
        else
            warning( ['Weights are < realmin from k about ' num2str(k) '.'] );
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the expansion of the orthonormal polynomial in the bulk without e^(x/2)
function p = pb(np, y, alpha, T)
z = y/4/np;
m2nxi = 2i*np*( sqrt(z).*sqrt(1 - z) - acos(sqrt(z) ) ); % = -2*n*xin
d = z - 1;
phi = 2*z - 1 + 2*sqrt(z)*sqrt(d);
if T == 1
    p = real( 1/z^(1/4)/(z-1)^(1/4)*...
        (exp(-m2nxi)*sqrt(phi)*(phi/z)^(alpha/2) + ...
        z^(-alpha)*exp(m2nxi)*1i/sqrt(phi)*(phi/z)^(-alpha/2) ) );
    return;
end
% Getting the higher order terms is hard-coded for speed and code length, but
% can be made to get arbitrary orders for general weight functions.
R1 = 0.0;
R2 = 0i;
if (alpha == 0)
    if ( T >= 5 )
        R1 = R1 + (-0.0001123610837959949*z^(-1) -1.490116119384768e-05*z^(-2) )/np^4;
        R2 = R2 + (-0.0002556506498360341i*z^(-1) -0.000452995300292969i*z^(-2) )/np^4;
        R1 = R1 + (+0.0001123610837959949*d^(-1) -4.159894009185921e-05*d^(-2) )/np^4;
        R2 = R2 + (+3.220671979488066e-05i*d^(-1) +1.006970189726202e-05i*d^(-2) )/np^4;
        R1 = R1 + (+0.0001014470072930732*d^(-3) -0.0001985441019505633*d^(-4) )/np^4;
        R2 = R2 + (-0.000150605189947433i*d^(-3) +0.004846322487411191i*d^(-4) )/np^4;
        R1 = R1 + (-0.0008181122595390668*d^(-5) -0.000793068006696034*d^(-6) )/np^4;
        R2 = R2 + (+0.02270678924434961i*d^(-5) +0.01903363216070482i*d^(-6) )/np^4;
    end
    if ( T >= 4 )
        R1 = R1 + (+0.0002585517035590279*z^(-1) +0.0005722045898437501*z^(-2) )/np^3;
        R2 = R2 + (-0.0009774102105034725i*z^(-1) +0.0005722045898437501i*z^(-2) )/np^3;
        R1 = R1 + (-0.0002585517035590275*d^(-1) -1.465597270447586e-05*d^(-2) )/np^3;
        R2 = R2 + (-0.000218577443817516i*d^(-1) -8.356541763117039e-05i*d^(-2) )/np^3;
        R1 = R1 + (+0.0003698466736593355*d^(-3) -0.002262821903935187*d^(-4) )/np^3;
        R2 = R2 + (+0.006447629575376164i*d^(-3) +0.02017682864342207i*d^(-4) )/np^3;
        R1 = R1 + (-0.008014160909770449*d^(-5) )/np^3;
        R2 = R2 + (+0.008014160909770449i*d^(-5) )/np^3;
    end
    if ( T >= 3 )
        R1 = R1 + (+0.0001627604166666672*z^(-1) )/np^2;
        R2 = R2 + (+0.004557291666666669i*z^(-1) )/np^2;
        R1 = R1 + (-0.0001627604166666668*d^(-1) -0.0007052951388888891*d^(-2) )/np^2;
        R2 = R2 + (-0.001085069444444443i*d^(-1) +0.01323784722222223i*d^(-2) )/np^2;
        R1 = R1 + (-0.001898871527777778*d^(-3) )/np^2;
        R2 = R2 + (+0.02278645833333334i*d^(-3) )/np^2;
    end
    R1 = R1 + (-0.015625*z^(-1) )/np^1;
    R2 = R2 + (-0.015625i*z^(-1) )/np^1;
    R1 = R1 + (+0.015625*d^(-1) -0.02604166666666667*d^(-2) )/np^1 + 1;
    R2 = R2 + (+0.05729166666666667i*d^(-1) +0.02604166666666667i*d^(-2) )/np^1;
elseif (alpha == 1)
    if ( T >= 5 )
        R1 = R1 + (-4.915484675654808e-05*z^(-1) +0.0003713369369506837*z^(-2) )/np^4;
        R2 = R2 + (-0.0003443859241626885i*z^(-1) +0.0002336502075195314i*z^(-2) )/np^4;
        R1 = R1 + (+4.915484675654708e-05*d^(-1) -6.56338876166932e-05*d^(-2) )/np^4;
        R2 = R2 + (-1.047197192785386e-05i*d^(-1) -0.0001160778626492966i*d^(-2) )/np^4;
        R1 = R1 + (+0.000152441503579723*d^(-3) -0.004777727892369407*d^(-4) )/np^4;
        R2 = R2 + (+0.004768120801007313i*d^(-3) +0.02510245091630599i*d^(-4) )/np^4;
        R1 = R1 + (-0.0228570547614078*d^(-5) -0.01982670016740085*d^(-6) )/np^4;
        R2 = R2 + (+0.0302200650972594i*d^(-5) +0.009516816080352408i*d^(-6) )/np^4;
    end
    if ( T >= 4 )
        R1 = R1 + (+0.0005976359049479174*z^(-1) -0.0008010864257812502*z^(-2) )/np^3;
        R2 = R2 + (+0.0009695688883463542i*z^(-1) -0.0002002716064453126i*z^(-2) )/np^3;
        R1 = R1 + (-0.0005976359049479167*d^(-1) +0.0008296636887538567*d^(-2) )/np^3;
        R2 = R2 + (-0.000400510246371044i*d^(-1) +0.006527672284915132i*d^(-2) )/np^3;
        R1 = R1 + (-0.007273111225646224*d^(-3) -0.01923398618344908*d^(-4) )/np^3;
        R2 = R2 + (+0.02336585433394821i*d^(-3) +0.01777258037049094i*d^(-4) )/np^3;
        R1 = R1 + (-0.008014160909770449*d^(-5) )/np^3;
        R2 = R2 + (+0.002003540227442612i*d^(-5) )/np^3;
    end
    if ( T >= 3 )
        R1 = R1 + (-0.00048828125*z^(-1) )/np^2;
        R2 = R2 + (-0.001953125000000002i*z^(-1) )/np^2;
        R1 = R1 + (+0.0004882812499999989*d^(-1) -0.01177300347222222*d^(-2) )/np^2;
        R2 = R2 + (+0.01323784722222223i*d^(-1) +0.02886284722222222i*d^(-2) )/np^2;
        R1 = R1 + (-0.02468532986111112*d^(-3) )/np^2;
        R2 = R2 + (+0.01139322916666667i*d^(-3) )/np^2;
    end
    R1 = R1 + (+0.046875*z^(-1) )/np^1;
    R2 = R2 + (+0.01171875i*z^(-1) )/np^1;
    R1 = R1 + (-0.046875*d^(-1) -0.02604166666666667*d^(-2) )/np^1 + 1;
    R2 = R2 + (+0.06119791666666667i*d^(-1) +0.006510416666666667i*d^(-2) )/np^1;
end
p = real( 1/z^(1/4)/(z-1+0i)^(1/4)*( (sqrt(phi)*R1 - 4^alpha*1i/sqrt(phi)*R2 ...
    )*(phi/z)^(alpha/2)*exp(-m2nxi) + (1i/sqrt(phi)*R1 + 4^alpha*sqrt(phi)*R2...
    )*(phi*z)^(-alpha/2)*exp(m2nxi) ) );

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the expansion of the orthonormal polynomial near zero without e^(x/2)
function p = pl(np, y, alpha, T)
% [FIXME] Ensure analytic continuation of all functions to avoid adding eps*1i.
z = (1+eps*1i)*y/4/np;
npb = 2*np*(pi/2 + sqrt(z).*sqrt(1 - z) - acos(sqrt(z) ) ); % = 2i n sqrtphitn

if T == 1
    p = real( sqrt(2*pi)*(-1)^np*sqrt(npb)/z^(1/4)/(1 - z)^(1/4)*...
        z^(-alpha/2)*(sin( (alpha + 1)/2*acos(2*z - 1) - pi*alpha/2)*...
        besselj(alpha,npb) + cos( (alpha + 1)/2*acos(2*z - 1) - ...
        pi*alpha/2)*(besselj(alpha-1,npb) - alpha/(npb)*...
        besselj(alpha, npb) ) ) );
    return
end
% Use the series expansion of R because it is faster and we use pl only very
% close to zero to have less calls to besselj.
R1 = 0i;
R2 = 0i;
if ( alpha == 0 )
    if ( T >= 5 )
        R1 = R1 + (+0.001501174033247706*z^4 +0.006473841918350666*z^3 )/np^4;
        R2 = R2 + (+0.9524842232974917i*z^4 +0.3550440568543063i*z^3 )/np^4;
        R1 = R1 + (+0.005174830777647138*z^2 +0.002518294998040369*z^1 )/np^4;
        R2 = R2 + (+0.1020198851574237i*z^2 +0.01827557013031548i*z^1 )/np^4;
        R1 = R1 + (+0.0006884162808641975 )/np^4;
        R2 = R2 + (+0.0009178883744855927i )/np^4;
    end
    if ( T >= 4 )
        R1 = R1 + (+0.8486152325952079*z^5 +0.4599905388515398*z^4 )/np^3;
        R2 = R2 + (+0.02362629629145646i*z^5 +0.07442783509439153i*z^4 )/np^3;
        R1 = R1 + (+0.2227888285411434*z^3 +0.0915386169900059*z^2 )/np^3;
        R2 = R2 + (+0.07472004547236026i*z^3 +0.05208755878894767i*z^2 )/np^3;
        R1 = R1 + (+0.02883322310405644*z^1 +0.005362654320987655 )/np^3;
        R2 = R2 + (+0.02605544532627866i*z^1 +0.008043981481481482i )/np^3;
    end
    if ( T >= 3 )
        R1 = R1 + (+0.03274424783742833*z^6 +0.02138984808385162*z^5 )/np^2;
        R2 = R2 + (-0.5345928880817182i*z^6 -0.3890960516718894i*z^5 )/np^2;
        R1 = R1 + (+0.01205177881103808*z^4 +0.004772192827748389*z^3 )/np^2;
        R2 = R2 + (-0.2664715497122905i*z^4 -0.1667621987066432i*z^3 )/np^2;
        R1 = R1 + (-0.0003747795414462069*z^2 -0.003240740740740739*z^1 )/np^2;
        R2 = R2 + (-0.09005731922398588i*z^2 -0.03657407407407406i*z^1 )/np^2;
        R1 = R1 + (-0.003472222222222222 )/np^2;
        R2 = R2 + (-0.006944444444444446i )/np^2;
    end
    R1 = R1 + (-0.2113083089153498*z^7 -0.1842728591546934*z^6 )/np^1;
    R2 = R2 + (+0.1378214216323702i*z^7 +0.1106475197021934i*z^6 )/np^1;
    R1 = R1 + (-0.1569653787325745*z^5 -0.129241088129977*z^4 )/np^1;
    R2 = R2 + (+0.08313723470337227i*z^5 +0.05509753620864732i*z^4 )/np^1;
    R1 = R1 + (-0.1008289241622575*z^3 -0.07116402116402117*z^2 )/np^1;
    R2 = R2 + (+0.02615520282186949i*z^3 -0.0044973544973545i*z^2 )/np^1;
    R1 = R1 + (-0.03888888888888889*z^1 -3.469446951953614e-18 )/np^1 + 1;
    R2 = R2 + (-0.0388888888888889i*z^1 -0.08333333333333333i )/np^1;
elseif ( alpha == 1 )
    
    if ( T >= 5 )
        R1 = R1 + (-0.9842194497752994*z^4 -0.3638886994207212*z^3 )/np^4;
        R2 = R2 + (-0.1111429035402386i*z^4 -0.06729307069816329i*z^3 )/np^4;
        R1 = R1 + (-0.1025708529296493*z^2 -0.01716940402704293*z^1 )/np^4;
        R2 = R2 + (-0.02748385600382132i*z^2 -0.005803066211052317i*z^1 )/np^4;
        R1 = R1 + (-0.0002294720936213962 )/np^4;
        R2 = R2 + (+0.0003351658950617355i )/np^4;
    end
    if ( T >= 4 )
        R1 = R1 + (+0.07190356645934227*z^5 -0.008994495608252214*z^4 )/np^3;
        R2 = R2 + (+0.2611881760337451i*z^5 +0.1415801037403221i*z^4 )/np^3;
        R1 = R1 + (-0.03246736545347658*z^3 -0.0269951499118166*z^2 )/np^3;
        R2 = R2 + (+0.06273458794292135i*z^3 +0.01863839285714289i*z^2 )/np^3;
        R1 = R1 + (-0.01297949735449737*z^1 -0.002681327160493828 )/np^3;
        R2 = R2 + (+0.001124338624338644i*z^1 -0.0004340277777777698i )/np^3;
    end
    if ( T >= 3 )
        R1 = R1 + (+0.5818842080253015*z^6 +0.422774248874778*z^5 )/np^2;
        R2 = R2 + (-0.1158436958397275i*z^6 -0.06655856841571131i*z^5 )/np^2;
        R1 = R1 + (+0.2885276040831597*z^4 +0.1792151675485009*z^3 )/np^2;
        R2 = R2 + (-0.02893518518518523i*z^4 -0.003174603174603194i*z^3 )/np^2;
        R1 = R1 + (+0.09497354497354499*z^2 +0.03611111111111111*z^1 )/np^2;
        R2 = R2 + (+0.01021825396825397i*z^2 +0.009722222222222219i*z^1 )/np^2;
        R1 = R1 + (+0.003472222222222224 )/np^2;
        R2 = R2 + (-0.01041666666666667i )/np^2;
    end
    R1 = R1 + (-0.147667867682724*z^7 -0.1203555135830268*z^6 )/np^1;
    R2 = R2 + (-0.01071718871208446i*z^7 -0.01703582895646388i*z^6 )/np^1;
    R1 = R1 + (-0.09264242400750337*z^5 -0.06428731762065096*z^4 )/np^1;
    R2 = R2 + (-0.0231811821335631i*z^5 -0.02897306397306398i*z^4 )/np^1;
    R1 = R1 + (-0.03481481481481481*z^3 -0.003174603174603177*z^2 )/np^1;
    R2 = R2 + (-0.03396825396825398i*z^3 -0.03690476190476191i*z^2 )/np^1;
    R1 = R1 + (+0.03333333333333333*z^1 +0.08333333333333333 )/np^1 + 1;
    R2 = R2 + (-0.03333333333333333i*z^1 -0.125i )/np^1;
end

p = real( sqrt(2*pi)*(-1)^np*sqrt(npb)/z^(1/4)/ ...
    (1 - z)^(1/4)*z^(-alpha/2)*( (sin( (alpha + 1)/2*acos(2*z - 1) - ...
    pi*alpha/2)*R1 -1i*sin( (alpha - 1)/2*acos(2*z - 1) - pi*alpha/2)*...
    R2*4^alpha)*besselj(alpha, npb) + (cos( (alpha + 1)/2*acos(2*z - ...
    1)- pi*alpha/2)*R1 - 1i*cos( (alpha - 1)/2*acos(2*z - 1) - pi*alpha/2)*...
    R2*4^alpha)*(besselj(alpha-1, npb) - ...
    alpha/npb*besselj(alpha, npb) ) ) );

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the expansion of the orthonormal polynomial near 4n without e^(x/2)
function p = pr(np, y, alpha, T)
z = y/4/np;
fn = (np*3*1i*( sqrt(z).*sqrt(1 - z) - acos(sqrt(z) ) ))^(2/3);
if T == 1
    p = real( 4*sqrt(pi)/z^(1/4)/d^(1/4)*z^(-alpha/2)* ...
        (cos( (alpha + 1)/2*acos(2*z - 1) )*fn^(1/4)*airy(0,fn) + ...
        -1i*sin( (alpha + 1)/2*acos(2*z - 1) )*fn^(-1/4)*airy(1,fn) ) );
    return
end
R1 = 0.0;
R2 = 0i;
d = z - 1;
if ( alpha == 0 )
    if ( T >= 5 )
        R1 = R1 + (-2.796022042244435e-05*d^4 +2.011295619003682e-05*d^3 )/np^4;
        R2 = R2 + (-0.002475382540803074i*d^4 +0.002022143985881445i*d^3 )/np^4;
        R1 = R1 + (-1.295611489560738e-05*d^2 +6.757898749300978e-06*d^1 )/np^4;
        R2 = R2 + (-0.001569093104321315i*d^2 +0.001116765488720164i*d^1 )/np^4;
        R1 = R1 + (-2.026867991649663e-06 )/np^4;
        R2 = R2 + (-0.000666664945504233i )/np^4;
    end
    if ( T >= 4 )
        R1 = R1 + (-0.003253020840537458*d^5 +0.002696059409107524*d^4 )/np^3;
        R2 = R2 + (-0.001993718633874486i*d^5 +0.001436399596567729i*d^4 )/np^3;
        R1 = R1 + (-0.002141460647635705*d^3 +0.001590456473343028*d^2 )/np^3;
        R2 = R2 + (-0.0008821041836667568i*d^3 +0.0003328986988216683i*d^2 )/np^3;
        R1 = R1 + (-0.001045363705958944*d^1 +0.0005110818194151536 )/np^3;
        R2 = R2 + (+0.000206871642466881i*d^1 -0.0007267403892403877i )/np^3;
    end
    if ( T >= 3 )
        R1 = R1 + (+7.18695260864394e-05*d^6 -6.692939666830364e-05*d^5 )/np^2;
        R2 = R2 + (+0.004480794599931599i*d^6 -0.004460948839100688i*d^5 )/np^2;
        R1 = R1 + (+6.109774830863335e-05*d^4 -5.407654074320793e-05*d^3 )/np^2;
        R2 = R2 + (+0.0044332342004771i*d^4 -0.004392654646940364i*d^3 )/np^2;
        R1 = R1 + (+4.5413316841889e-05*d^2 -3.439153439153486e-05*d^1 )/np^2;
        R2 = R2 + (+0.004329315657887089i*d^2 -0.004221019721019722i*d^1 )/np^2;
        R1 = R1 + (+1.984126984127046e-05 )/np^2;
        R2 = R2 + (+0.004007936507936511i )/np^2;
    end
    R1 = R1 + (+0.01272750016968636*d^7 -0.01255318649195314*d^6 )/np^1;
    R2 = R2 + (+0.01249169963638477i*d^7 -0.01227111943654102i*d^6 )/np^1;
    R1 = R1 + (+0.01234301148871137*d^5 -0.01208266968103703*d^4 )/np^1;
    R2 = R2 + (+0.0119973041110328i*d^5 -0.01164522908686174i*d^4 )/np^1;
    R1 = R1 + (+0.01174829614829615*d^3 -0.01129622758194187*d^2 )/np^1;
    R2 = R2 + (+0.01117002997002997i*d^3 -0.01048155019583591i*d^2 )/np^1;
    R1 = R1 + (+0.01063492063492063*d^1 -0.009523809523809525 )/np^1 + 1;
    R2 = R2 + (+0.009365079365079364i*d^1 -0.007142857142857144i )/np^1;
elseif ( alpha == 1 )
    if ( T >= 5 )
        R1 = R1 + (+0.001882808665107142*d^4 -0.001505607402838352*d^3 )/np^4;
        R2 = R2 + (+0.0007698736891847038i*d^4 -0.00053700578303835i*d^3 )/np^4;
        R1 = R1 + (+0.001127942655094661*d^2 -0.0007500605012473177*d^1 )/np^4;
        R2 = R2 + (+0.0003047918362998525i*d^2 -7.445206877896294e-05i*d^1 )/np^4;
        R1 = R1 + (+0.0003729230795475838 )/np^4;
        R2 = R2 + (-0.000150200751078725i )/np^4;
    end
    if ( T >= 4 )
        R1 = R1 + (+0.004744560105086139*d^5 -0.003932667625282309*d^4 )/np^3;
        R2 = R2 + (+0.0003596418409404624i*d^5 -0.0001584559939677844i*d^4 )/np^3;
        R1 = R1 + (+0.003122487606361571*d^3 -0.002315841844741004*d^2 )/np^3;
        R2 = R2 + (-4.037681508993723e-05i*d^3 +0.0002340840897612615i*d^2 )/np^3;
        R1 = R1 + (+0.001516718193622956*d^1 -0.000734979989146656 )/np^3;
        R2 = R2 + (-0.0004153882471144348i*d^1 +0.0005610623173123191i )/np^3;
    end
    if ( T >= 3 )
        R1 = R1 + (-0.0006821631748549443*d^6 +0.0006790103374176761*d^5 )/np^2;
        R2 = R2 + (-0.001811910501552051i*d^6 +0.001788299650036054i*d^5 )/np^2;
        R1 = R1 + (-0.0006715318212236979*d^4 +0.0006562001666763589*d^3 )/np^2;
        R2 = R2 + (-0.00175433095900883i*d^4 +0.001701906189049048i*d^3 )/np^2;
        R1 = R1 + (-0.0006256188256188247*d^2 +0.0005622895622895622*d^1 )/np^2;
        R2 = R2 + (-0.001612467162467158i*d^2 +0.001434463684463688i*d^1 )/np^2;
        R1 = R1 + (-0.0004166666666666676 )/np^2;
        R2 = R2 + (-0.0009722222222222267i )/np^2;
    end
    R1 = R1 + (-0.04460324252123175*d^7 +0.04443907573247043*d^6 )/np^1;
    R2 = R2 + (-0.01087128029713879i*d^7 +0.01076543684375636i*d^6 )/np^1;
    R1 = R1 + (-0.04423440188249792*d^5 +0.04396981497716192*d^4 )/np^1;
    R2 = R2 + (-0.01062009860203378i*d^5 +0.0104078944789149i*d^4 )/np^1;
    R1 = R1 + (-0.04361026909598338*d^3 +0.04308472479901051*d^2 )/np^1;
    R2 = R2 + (-0.01006952095523524i*d^3 +0.009451041022469598i*d^2 )/np^1;
    R1 = R1 + (-0.04222222222222222*d^1 +0.04047619047619047 )/np^1 + 1;
    R2 = R2 + (-0.008015873015873025i*d^1 -0.03928571428571431i )/np^1;
end

p = real( 4*sqrt(pi)/z^(1/4)/d^(1/4)*z^(-alpha/2)* ...
    ( (R1*cos( (alpha + 1)/2*acos(2*z - 1) ) -1i*cos( (alpha - 1)/2* ...
    acos(2*z - 1) )*R2*4^alpha )*fn^(1/4)*airy(0,fn) + ...
    (-1i*sin( (alpha + 1)/2*acos(2*z - 1) )*R1 -sin( (alpha - 1)/2*...
    acos(2*z - 1) )*R2*4^alpha)*fn^(-1/4)*airy(1,fn) ) );

end
