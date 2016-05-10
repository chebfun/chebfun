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
%   By default LAGPTS uses 'GW' when N < 128.
%
% References:
%   [1] G. H. Golub and J. A. Welsch, "Calculation of Gauss quadrature rules",
%       Math. Comp. 23:221-230, 1969,
%   [2] A. Glaser, X. Liu and V. Rokhlin, "A fast algorithm for the calculation
%       of the roots of special functions", SIAM Journal on Scientific 
%       Computing", 29(4):1420-1438:, 2007.
%
% See also CHEBPTS, LEGPTS, HERMPTS, and JACPTS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
%
% 'GW' by Nick Trefethen, March 2009 - algorithm adapted from [1].
% 'FAST' by Nick Hale, March 2010 - algorithm adapted from [2].

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

% Check the inputs
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
    if ( ~any(strcmpi(method, {'default', 'GW', 'fast'})) )
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

% decide to use GW or FAST
if ( (n < 128 || strcmpi(method,'GW')) && ( ~strcmpi(method,'fast'))  )
    % GW, see [1]
    
    alpha = 2*(1:n)-1;  beta = 1:n-1;     % 3-term recurrence coeffs
    T = diag(beta,1) + diag(alpha) + diag(beta,-1);  % Jacobi matrix
    [V, D] = eig(T);                      % eigenvalue decomposition
    [x, indx] = sort(diag(D));            % Laguerre points
    w = V(1,indx).^2;                     % Quadrature weights
    v = sqrt(x).*abs(V(1,indx)).';        % Barycentric weights
    v = v./max(v); 
    v(2:2:n) = -v(2:2:n);
    
else
    % Fast, see [2]
    [x, ders] = alg0_Lag(n);              % Nodes and L_n'(x)
    w = exp(-x)./(x.*ders.^2); w = w';    % Quadrature weights
    v = exp(-x/2)./ders;                  % Barycentric weights
    v = -v./max(abs(v));
    
end
w = (1/sum(w))*w;                         % Normalise so that sum(w) = 1

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

[ignored, d1] = eval_Lag(x1, n);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Evauate Laguerre polynomial via recurrence
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
