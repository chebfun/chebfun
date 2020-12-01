function [x, w] = pswfpts(N, c, dom, quadtype)
%PSWFPTS   Quadrature nodes and weights from PSWF roots.
% X = PSWFPTS(N, C) returns the N roots of the Nth PSWF with bandwidth C.
% (See HELP PSWF for the definition.) 
%   
% [X, W] = PSWFPTS(N, C) also returns the weights for the interpolatory PSWF
% quadrature rule with the nodes X.
%
% [X, W] = PSWFPTS(N, C, DOM) scales the nodes and weights to the interval DOM,
% which should be a finite 2-vector.
%
% [X, W] = PSWFPTS(N, C, DOM, 'GGQ')  or PSWFPTS(N, C, 'GGQ') returns
% rather the nodes and weights corresponding to the N-point generalised
% Gauss quadrature rule, which is exact for PSWFs with bandwidth C of order
% up to 2N-1.
%
% Example:
%  f = pswf(8,pi); sum(f)
%  [x,w] = pswfpts(5,pi,[-1,1],'GGQ'); w*f(x)
%
% See also PSWF, LEGPTS.

% Copyright 2020 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer note: The approach used is to compute the Legendre coefficients
% of the degree N PSWF using the CHEBFUN.PSWF code and then find the
% roots with the Legendre analogue of the colleague matrix [1]. This is OK
% for small N and C, but for larger values more advanced techniques should
% be used; for instance as described in [2].
%
% For GGQ we use the ideas from [3]. See further comments in the subroutine.
%
% [1] R. M. Corless and G. Litt, Generalized companion matrices for polynomials
% not expressed in monomial Bases (Unpublished note)
% [2] A. Glaser, X. Liu, and V. Rokhlin, A fast algorithm for the roots of
% special functions, SISC, 29, 4, pp. 1420-1438, 2007.
% [3] J. Ma, V. Rokhlin, and S. Wandzura, "Generalised Gaussian Quadrature
% rules For systems of arbitrary functions", SINUM, 1996.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Defaults:
if ( nargin < 3 )
    dom = [-1 1];
    quadtype = 'roots';
end
if ( nargin == 3 )
    if ( isnumeric(dom) )   
        quadtype = 'roots';
    else
        quadtype = dom;
        dom = [-1 1];
    end
end

% Parse inputs:
assert( nargin >= 2, 'PSWFPTS requires at least two input arguments.')
assert( (numel(N)==1) && (round(N)==N) && (N>=0) , ...
    'Input N must be a non-negative integer.');
assert( (numel(c)==1) && (c>0) , ...
    'Input C must be a non-negative scalar.');
assert( numel(dom)==2 && all(isfinite(dom)) , ...
    'Domain must be a finite two-vector.');

% Trivial case:
if ( N == 0 )
    x = []; w = [];
    return
end

% Fork based on desired quadratyure type (PSWF roots or GGQ):
if ( strcmpi(quadtype, 'ggq') )
    [x,w] = pswfggq(N, c);
else
    [x,w] = pswfrootsquad(N, c);
end

% Enforce symmetry:
x = mean([x -flipud(x)],2);
w = mean([w; fliplr(w)]);
    
% Scale if required:
if ( nargin >= 3 )
    x = scaleNodes(x, dom);
    w = scaleWeights(w, dom);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,w] = pswfrootsquad(N, c)
%ROOTSQUAD   Quadrature rule with nodes at PSWF roots.
% [X,W] = ROOTSQUAD(N,C) returns the N roots the Nth PSWF with bandwidth C
% and the weights for the associated interpolatory quadrature rule that
% integrates PSWFs with bandwidth C and degree 0...N-1 exactly.

% Obtain Legendre coeffs of PSWFs of up to degree N:
V = pswf(0:N, c, [-1 1], 'coeffs');

% Compute the roots of P_{N+1}:
x = legroots(V(:,end));

if ( nargout > 1 )
    % Integrate the PSWFs (use Legendre orthogonality):
    S = 2*V(1,1:N); 
    % Construct the Legendre-Vandermonde matrix:
    P = legpolyval(V(:,1:N), x);
    % Compute the weights:
    w = S/P;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,w] = pswfggq(N, c, x)
%Generalised (PSWF) Gauss quadrature.
% [X, W] = ggq(N, C) returns the N quadrature nodes and weights for
% the generalised Gauss quadrature method which exactly integrates PSWFs
% of orders 0 to 2N-1 with bandwidth parameter C.

% See [1] Ma, Rokhlin, Wandzura, "Generalised Gaussian Quadrature Rules For
% Systems of Arbitrary Functions", SINUM 1996.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer note: The modified Newton method is very sensitive to the
% intial guess. Some remedies are suggested in [1], but here we just choose
% what we hope is a good guess and cross our fingers. It seems to work well
% for smallish N and C.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( nargin == 2 )
    % Initial guess:
    % x = pswfpts(N, 2*c); % PSWF roots
    x = legpts(N); % Gauss-Legendre nodes
    % Use KTE map to improve chances of Newton convergence!
    a = .5; x = asin(a*x)./asin(a); % Hello, my old friend!
end

% Continuation for c > N. (Hacky! should be improved.)
if ( nargin == 2 && c > N )
    for cc = [N:2:c c]
        [x, w] = pswfggq(N, cc, x);
    end
    return
end

% Get Legendre coeffs of the first 2N PSWFs:
V = pswf(0:2*N-1, c, [-1 1], 'coeffs'); 
S = 2*V(1,:); % Integral of each column (using orthogonality of Legendre)
% Differentiate the PSWFs:
M = size(V,1);
C = ultraS.convertmat(M-1, 0.5, 0.5);
Vp = C\V(2:end,:);

% Modified Newton iteration to find x (see [1]):
A = zeros(2*N,2*N);
for k = 1:10 % Quadratic convergence expected, so 10 iterations should be OK.
    % Construct Legendre-Hermite-Vandermonde matrix:  
    A(1:2:end,:) = legpolyval(V,x);
    A(2:2:end,:) = legpolyval(Vp,x);
    % Invert and integrate (can be reduced to one solve):
    w = S/A;
    % Modified Newton update.
    dx = w(2:2:end)./w(1:2:end);
    x = x + dx.';
    if ( norm(dx, inf) < 1e-12 ), break, end % Escape clause
end

if ( any(isnan(x)) || any(abs(x) > 1) )
    error('CHEBFUN:pswfpts:iterationfailure', ...
        'Newton iteration failed to converge.');
elseif ( (norm(dx, inf) > 1e-10) )
    warning('CHEBFUN:pswfpts:possibleiterationfailure', ...
        'Newton iteration may have failed to converge.');
end

% Compute the weights:
w = S/legpolyval(V,x);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = legroots(v)
%LEGROOTS    Roots of a Legendre series
% X = LEGROOTS(V) returns the real-valued roots in [-1, 1] of the
% polynomial whose Legendre series coefficients are given by the vector V.
% Note that the coefficients in V should be given in DESCENDING order,
% i.e., P(x) = C(1)*P_0(x) + C(2)*P_1(x) + ... + C(N)*P_{N-1}(x)
%
% See also CHEBTECH/ROOTS.

% The roots are found by solving a 'companion matrix' eigenvalue problem.
% See RM Corless and G Litt, "Generalized Companion Matrices for Polynomials
% not expressed in Monomial Bases" (Unpublished note)

% Remove trailing zeros:
v = v(1:find(v, 1, 'last'));

% Legendre recursion coefficients:
n = length(v)-1; N = 0:n; 
a = (N+1)./(2*N+1); g = N./(2*N+1);

% Construct matrices:
C = diag(a(1:n-1), 1) + diag(g(2:n), -1);
C(end,:) = -v(1:n)';
D = -v(n-1)+(1-1/n)*v(n+1);
C(end,end-1) = D;
B = eye(n); B(n,n) = (2-1/n)*v(n+1);

% Generalised eigenvalue problem:
x = eig(C, B);

% Remove spurious values:
x(abs(imag(x)) > 1e-13) = [];
x(abs(x) > 1) = [];
x = sort(x);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = legpolyval(c, x)
%LEGPOLYVAL   Evaluate a Legendre polynomial.
% Y = LEGPOLYVAL(C, X), when C is a column vector of length N+1 whose
% elements are the Legendre coefficients of a polynomial, is the value of the
% polynomial evaluated at X, i.e.,
% 
%   Y = C(1)*P_N(X) + C(2)*T_{N-1}(X) + ... + C(N)*T_1(X) + C(N+1)*I
% 
% If C is an (N+1) x M matrix, then LEGPOLYVAL interprets each of the columns
% of C as coefficients of a degree N polynomial and evaluates the M Legendre
% expansions
% 
%   Y_m = C(1,m)*T_N(X) + ... + C(N,m)*T_1(X) + C(N+1,m)*P_0(X), 1 <= m <= M,
% 
% returning the results as columns of a matrix Y = [Y_1 ... Y_M].
%
% X must be a scalar or a column vector.

% Modified Clenshaw scheme:
x = repmat(x(:), 1, size(c, 2));
bk1 = zeros(size(x, 1), size(c, 2)); 
bk2 = bk1;
e = ones(size(x, 1), 1);
n = size(c,1)-1;
for k = n:-1:1
    bk = e*c(k+1,:) + (2*k+1)/(k+1)*x.*bk1 - (k+1)/(k+2)*bk2;
    bk2 = bk1; 
    bk1 = bk;
end
y = e*c(1,:) + x.*bk1 - .5*bk2;
end

% function y = legpolyval(c, x)
% %LEGPOLYVAL   Evaluate a Legendre polynomial.
% % Y = LEGPOLYVAL(C, X), when C is a column vector of length N+1 whose
% % elements are the Legendre coefficients of a polynomial, is the value of the
% % polynomial evaluated at X, i.e.,
% % 
% %   Y = C(1)*P_N(X) + C(2)*T_{N-1}(X) + ... + C(N)*T_1(X) + C(N+1)*I
% % 
% % If C is an (N+1) x M matrix, then LEGPOLYVAL interprets each of the columns
% % of C as coefficients of a degree N polynomial and evaluates the M Legendre
% % expansions
% % 
% %   Y_m = C(1,m)*T_N(X) + ... + C(N,m)*T_1(X) + C(N+1,m)*P_0(X), 1 <= m <= M,
% % 
% % returning the results as columns of a matrix Y = [Y_1 ... Y_M].
% %
% % X must be a scalar or a column vector.
%
% Construct Legendre-Vandermonde matrix and multiply. (This is faster, but
% possibly less accurate tan using the modified Clenshaw scheme above.)
% x = x(:);
% n = size(c,1)-1;
% P = zeros(length(x),n+1);
% P(:,1:2) = [1+0*x, x];
% for k = 1:n-1
%     P(:,k+2) = ((2*k+1)*x.*P(:,k+1) - k*P(:,k))/(k+1);
% end
% y = P*c;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = scaleNodes(x, dom)
%SCALENODES   Scale the nodes X from [-1,1] to DOM.
% TODO: Deal with unbounded domains
if ( dom(1) == -1 && dom(2) == 1 )
    % Nodes are already on [-1, 1];
    y = x;
    return
end
% Scale the nodes:
y = dom(2)*(x + 1)/2 + dom(1)*(1 - x)/2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w = scaleWeights(w, dom)
%SCALEWEIGHTS   Scale the weights W from [-1,1] to DOM.
% TODO: Deal with unbounded domains
if ( dom(1) == -1 && dom(2) == 1 )
    % Nodes are already on [-1, 1];
    return
end
% Scale the weights:
w = (diff(dom)/2)*w;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
