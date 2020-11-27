function [x, w] = pswfpts(N, c, dom, quadtype)
%PSWFPTS   Quadrature nodes and weights from PSWF roots.
% X = PSWFPTS(N, C) returns the N roots of the (N+1)st prolate spheroidal
% wave function with bandwidth C.
%   
% [X, W] = PSWFPTS(N, C) returns also the weights for the interpolatory PSWF
% quadrature rule with the nodes X.
%
% [X, W] = PSWFPTS(N, C, DOM) scales the nodes and weights to the interval DOM,
% which should be a finite two-vector.
%
% [X, W] = PSWFPTS(N, C, DOM, 'GGQ') returns rather the nodes and weights
% corresponding to the N-point generalised Gauss quadrature rule, which is
% exact for PSWFs with bandwidth C of order up to 2N.
%
% Example:
%
% f = pswf(9,pi); sum(f)
% [x,w] = pswfpts(5,pi,[-1,1],'GGQ'); w*f(x)
%
% See also PSWF, LEGPTS.

% Copyright 2020 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer note: The approach used is to compute the Legendre coefficients
% of the degree N+1 PSWF using the CHEBFUN.PSWF code and then find the
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
%
% Additional comment: Symmetry is not enforced in either the nodes or weights.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( nargin == 1 )
    c = N;
end
if ( nargin < 4 )
    quadtype = 'roots';
end    

if ( strcmpi(quadtype, 'ggq') )
    [x,w] = ggq(N, c);
else
    [x,w] = rootsquad(N, c);
end
    
if ( nargin >= 3 )
    x = scaleNodes(x, dom);
    w = scaleWeights(w, dom);
end

end

function [x,w] = rootsquad(N, c)
%ROOTSQUAD   Quadrature rule with nodes at PSWF roots.
% [X,W] = ROOTSQUAD(N,C) returns the N roots the (N+1)st PSWF with
% bandwidth C and the weights for the associated interpolatory quadrature
% rule.

% Obtain Legendre coeffs of PSWFs of up to degree N+1:
V = pswf(1:N+1, c, [-1 1], 'coeffs');

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

function [x,w] = ggq(N, c)
%Generalised (PSWF) Gauss quadrature.
% [X, W] = ggq(N, C) returns the N quadrature nodes and weights for
% the generalised Gauss quadrature method which exactly integrates PSWFS
% of orders 1 to 2N with bandwidth parameter C.

% See [1] Ma, Rokhlin, Wandzura, "Generalised Gaussian Quadrature Rules For
% Systems of Arbitrary Functions", SINUM 1996.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer note: The modified Newton method is very sensitive to the
% intial guess. Some remedies are suggested in [1], but here we just choose
% what we hope is a good guess and cross our fingers. It seems to work well
% for smallish N and C.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial guess:
% x = pswfpts(N, c); % PSWF roots
x = legpts(N); % Gauss-Legendre nodes
% Use KTE map to improve chances of Newton convergence!
a = .5; x = asin(a*x)./asin(a); % Hello, my old friend!

% Get Legendre coeffs of the first 2N PSWFs:
V = pswf(1:2*N, c, [-1 1], 'coeffs'); S = 2*V(1,:);
% Differentiate the PSWFs:
M = size(V,1);
C = ultraS.convertmat(M-1, 0.5, 0.5);
D = ultraS.diffmat(M-1, 0);
Vp = C\(D*V(2:end,:));

% Modified Newton iteration to find x (see [1]):
A = zeros(2*N,2*N);
for k = 1:10
    % Construct Legendre-Hermite-Vandermonde matrix:
    A(1:2:end,:) = legpolyval(V,x);
    A(2:2:end,:) = legpolyval(Vp,x);
    w = S/A;
    dx = w(2:2:end)./w(1:2:end);
    x = x + dx.';
    if ( norm(dx, inf) < 1e-10 )
        break
    end
end

if ( (norm(dx) > 1e-10) || any(isnan(x)) )
    warning('CHEBFUN:pswfpts:iterationfailure', ...
        'Newton iteration may have failed to converge.');
end

% Compute the weights:
w = S/legpolyval(V,x);

% % Test by integating for k = 1..2n (TODO: Remove.)
% P = pswf(1:2*N, c, [-1 1]);
% err = norm(w*P(x) - integral(P), inf)

end


function x = legroots(v)
%LEGROOTS    Roots of a Legendre series

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

function y = legpolyval(c, x)
% Clenshaw scheme for array-valued functions.
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
% x = x(:);
% n = size(c,1)-1;
% P = zeros(length(x),n+1);
% P(:,1:2) = [1+0*x, x];
% for k = 1:n-1
%     P(:,k+2) = ((2*k+1)*x.*P(:,k+1) - k*P(:,k))/(k+1);
% end
% y = P*c;
% end

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
