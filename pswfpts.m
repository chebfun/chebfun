function [x, w] = pswfpts(N, c, dom)
%PSWFPTS   Quadrature nodes and weights from PSWF roots
% X = PSWFPTS(N, C) returns the N roots the (N+1)st PSWF with bandwidth C.
%   
% [X,W] = PSWF(N, C) returns also the weights for the interpolatory PSWF
% quadrature rule with the nodes X.
%
% See also PSWF, LEGPTS.

% Copyright 2020 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer note: The approach used is to compute the Legendre coefficients
% of the degree N+1 PSWF using the CHEBFUN.PSWF code and then find the
% roots with the Legendre analogue of the Colleague matrix [1]. This is OK
% for small N and C, but for larger values more advanced techniques should
% be used; for instance as described in [2].
%
% [1] RM Corless and G Litt, Generalized Companion Matrices for Polynomials
% not expressed in Monomial Bases (Unpublished note)
% [2] A Glaser, X Liu, and V Rokhlin, A fast algorithm for the roots of
% special functions, SISC, 29, 4, p1420-1438, 2007.
%
% Additional comment: Symmetry is not enforced in either the nodes or weights.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( nargin == 1 )
    c = N;
end

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

if ( nargin == 3 )
    x = scaleNodes(x, dom);
    if ( nargout > 1 )
        w = scaleWeights(w);
    end
end

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

function y = scaleNodes(x, dom)
%SCALENODES   Scale the Chebyshev nodes X from [-1,1] to DOM.
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
%SCALEWEIGHTS   Scale the Chebyshev weights W from [-1,1] to DOM.
% TODO: Deal with unbounded domains
if ( dom(1) == -1 && dom(2) == 1 )
    % Nodes are already on [-1, 1];
    return
end
% Scale the weights:
w = (diff(dom)/2)*w;
end