function [x w v] = hermpts(n,varargin)
%HERMPTS  Hermite points and Gauss-Hermite Quadrature Weights.
%  HERMPTS(N) returns N Hermite points X in (-inf,inf). By default these are
%  roots of the 'physicist'-type Hermite polynomials, which are orthogonal
%  with respect to the weight exp(-x.^2).
%
%  HERMPTS(N,'PROB') normalises instead by the probablist's definition
%  (with weight exp(-x.^2/2)), which gives rise to monomials.
%
%  [X,W] = HERMPTS(N) returns also a row vector W of weights for Gauss-Hermite
%  quadrature. [X,W,V] = HERMPTS(N) returns in addition a column vector V
%  of the barycentric weights corresponding to X.
%
%  [X,W] = HERMPTS(N,METHOD) where METHOD is one of 'GW' or 'FAST' allows 
%  the user to select which method is used. 'GW' will use the traditional 
%  Golub-Welsch eigenvalue method and is best when N is small. 'FAST' will 
%  uses Glaser-Liu-Rokhlin fast algorithm which is much faster for large N.
%  By default HERMPTS uses 'GW' when N < 128.
%
%  See also chebpts, legpts, lagpts, and jacpts.

%  Copyright 2011 by The University of Oxford and The Chebfun Developers. 
%  See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

%  'GW' by Nick Trefethen, March 2009 - algorithm adapted from [1].
%  'FAST' by Nick Hale, March 2010 - algorithm adapted from [2].
%
%  References:
%   [1] G. H. Golub and J. A. Welsch, "Calculation of Gauss quadrature
%       rules", Math. Comp. 23:221-230, 1969,
%   [2] A. Glaser, X. Liu and V. Rokhlin, "A fast algorithm for the
%       calculation of the roots of special functions", SIAM Journal
%       on Scientific Computing", 29(4):1420-1438:, 2007.

% Defaults
method = 'default';
type = 'phys';

if n < 0
    error('CHEBFUN:hermpts:n','First input should be a positive number');
end

% Return empty vector if n == 0
if n == 0
    x = []; w = []; v = []; return
end

% Check the inputs
while ~isempty(varargin)
    s = varargin{1}; varargin(1) = [];
    if strcmpi(s,'GW'), method = 'GW';
    elseif strcmpi(s,'fast'), method = 'fast';   
    elseif strncmpi(s,'phys',3), type = 'phys'; 
    elseif strncmpi(s,'prob',3), type = 'prob'; 
    else error('CHEBFUN:hermpts:input',['Unrecognised input string.', s]); 
    end
end


% Decide to use GW or FAST
if n == 1
    x = 0; w = sqrt(pi); v = 1;           % n = 1 case is trivial
elseif (n < 128 || strcmpi(method,'GW')) && ~strcmpi(method,'fast') % GW, see [1]
    beta = sqrt(.5*(1:n-1));              % 3-term recurrence coeffs
    T = diag(beta,1) + diag(beta,-1);     % Jacobi matrix
    [V,D] = eig(T);                       % Eigenvalue decomposition
    [x,indx] = sort(diag(D));             % Hermite points
    w = sqrt(pi)*V(1,indx).^2;            % weights
    v = abs(V(1,indx)).';                 % Barycentric weights
    v = v./max(v); v(2:2:n) = -v(2:2:n);
    
    % Enforce symmetry
    ii = 1:floor(n/2);  x = x(ii);  w = w(ii);
    vmid = v(floor(n/2)+1); v = v(ii);
    if mod(n,2)
        x = [x ; 0 ; -flipud(x)];   w = [w  sqrt(pi)-sum(2*w) fliplr(w)];
        v = [v ; vmid ; v(end:-1:1)];
    else
        x = [x ; -flipud(x)];       w = [w fliplr(w)];
        v = [v ; -v(end:-1:1)];
    end
else                                                            % Fast, see [2]
    [x ders] = alg0_Herm(n);              % Nodes and H_n'(x)
    w = (2*exp(-x.^2)./ders.^2)';         % Quadrature weights
    v = exp(-x.^2/2)./ders;               % Barycentric weights
    v = v./max(abs(v));
    if ~mod(n,2), ii = (n/2+1):n; v(ii) = -v(ii); end
end
w = (sqrt(pi)/sum(w))*w;                  % Normalise so that sum(w) = sqrt(pi)

if strcmpi(type,'prob')
    x = x*sqrt(2);
    w = w*sqrt(2);
end

% -------------------- Routines for FAST algorithm ------------------------

function [roots ders] = alg0_Herm(n) % driver for 'Fast'.
% Compute coefficients of P_m(0), Pm'(0), m = 0,..,N
Hm2 = 0; Hm1 = pi^(-1/4); Hpm2 = 0; Hpm1 = 0;
for k = 0:n-1
    H = -sqrt(k/(k+1))*Hm2;
    Hp = sqrt(2/(k+1))*Hm1-sqrt(k/(k+1))*Hpm2;
    Hm2 = Hm1; Hm1 = H; Hpm2 = Hpm1; Hpm1 = Hp;
end

roots = zeros(n,1); ders = zeros(n,1);                      % allocate storage
if mod(n,2),roots((n-1)/2) = 0; ders((n+1)/2) = Hp;         % zero is a root
else [roots(n/2+1) ders(n/2+1)] = alg2_Herm(H,n); end        % find first root

[roots ders] = alg1_Herm(roots,ders); % compute roots and derivatives

% -------------------------------------------------------------------------

function [roots ders] = alg1_Herm(roots,ders) % main algorithm for 'Fast'
n = length(roots);
if mod(n,2), N = (n-1)/2; s = 1;
else N = n/2; s = 0; end

m = 30; % number of terms in Taylor expansion
hh1 = ones(m+1,1); u = zeros(1,m+1); up = zeros(1,m+1);
for j = N+1:n-1
    x = roots(j); % previous root
    
    % initial approx
    h = rk2_Herm(pi/2,-pi/2,x,n) - x;

    % scaling
    M = 1/h;
    
    % recurrence relation for Hermite polynomials
    c1 = -(2*n+1-x^2)/M^2; c2 = 2*x./M^3; c3 = 1./M^4;
    u(1) = 0; u(2) = ders(j)/M; 
    u(3) = .5*c1*u(1);
    u(4) = (c1*u(2) + c2*u(1))/6;
    up(1) = u(2); up(2) = 2*u(3)*M; up(3) = 3*u(4)*M; up(m+1) = 0;
    
    for k = 2:m-2
        u(k+3) = (c1*u(k+1) + c2*u(k) + c3*u(k-1))/((k+1)*(k+2));
        up(k+2) = (k+2)*u(k+3)*M;
    end
  
    % flip for more accuracy in inner product calculation
    u = u(m+1:-1:1);       up = up(m+1:-1:1);
    
    % Newton iteration
    hh = hh1; hh(end) = M;    step = inf;  l = 0; 
    while (abs(step) > eps) && (l < 10)
        l = l + 1;
        step = (u*hh)/(up*hh);
        h = h - step;
        hh = [M;cumprod(M*h+zeros(m,1))]; % powers of h (This is the fastest way!)
        hh = hh(end:-1:1); % flip for more accuracy in inner product calculation
    end
    
    % update
    roots(j+1) = x + h;
    ders(j+1) = up*hh;
end

% nodes are symmetric
roots(1:N+s) = -roots(n:-1:N+1);
ders(1:N+s) = ders(n:-1:N+1);

function [x1 d1] = alg2_Herm(Hn0,n) % find the first root (note H_n'(0)==0)

% advance ODE via Runge-Kutta for initial approx
x1 = rk2_Herm(0,-pi/2,0,n);

m = 30; % number of terms in Taylor expansion

% scaling
M = 1/x1;
% c = log10(n);
% M = 1./x1.^(1-1.25/(c));

% recurrence relation for Legendre polynomials
u = zeros(1,m+1); up = zeros(1,m+1);
u(1) = Hn0; u(3) = -.5*(2*n+1)*u(1)/M^2;
up(1) = 0; up(2) = 2*u(3)*M;
for k = 2:2:m-2
    u(k+3) = (-(2*n+1)*u(k+1)/M^2 + u(k-1)/M^4)/((k+1)*(k+2));
    up(k+2) = (k+2)*u(k+3)*M;
end

% flip for more accuracy in inner product calculation
u = u(m+1:-1:1);
up = up(m+1:-1:1);

x1k = [M;cumprod(M*x1+zeros(m,1))];

step = inf; l = 0;
% Newton iteration
while (abs(step) > eps) && (l < 10)
    l = l + 1;
    step = (u*x1k)/(up*x1k);
    x1 = x1 - step;
    x1k = [1;cumprod(M*x1+zeros(m,1))]; % powers of h (This is the fastest way!)
    x1k = x1k(end:-1:1);
end
% Update derivative
d1 = up*x1k;

% -------------------------------------------------------------------------

function x = rk2_Herm(t,tn,x,n) % Runge-Kutta for Hermite Equation
m = 10; h = (tn-t)/m;
for j = 1:m
    k1 = -h/(sqrt(2*n+1-x^2)-.5*x*sin(2*t)/(2*n+1-x^2));
    t = t+h;
    k2 = -h/(sqrt(2*n+1-(x+k1)^2)-.5*x*sin(2*t)/(2*n+1-(x+k1)^2));
    x = x+.5*(k1+k2);
end

