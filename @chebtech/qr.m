function [f, R, E] = qr(f, outputFlag, methodFlag)
%QR   QR factorisation of an array-valued CHEBTECH.
%   [Q, R] = QR(F) returns a QR factorisation of F such that F = Q*R, where the
%   CHEBTECH Q is orthogonal (with respect to the continuous L^2 norm on [-1,1])
%   and of the same size as F and R is an m x m upper-triangular matrix when F
%   has m columns.
%
%   [Q, R, E] = QR(F) produces unitary Q, upper-triangular R, and a permutation
%   matrix E so that F*E = Q*R. The column permutation E is chosen to reduce
%   fill-in in R.
%
%   [Q, R, E] = QR(F, 'vector') returns the permutation information as a vector
%   instead of a matrix.  That is, E is a row vector such that F(:,E) = Q*R.
%   Similarly, [Q, R, E] = QR(F, 'matrix') returns a permutation matrix E. This
%   is the default behavior.
%
%   QR(F, 'vector', METHOD) or QR(F, 'vector', METHOD) specifies which method
%   to use in computing the QR factorisation. METHOD = 'built-in' will for a
%   weighted Legendre-Vandermonde matrix and orthogonalise this with the
%   standard Matlab QR algorithm. METHOD = 'householder' uses the technique
%   described in [1]. METHOD = 'householder' is the default option.
%
%   [1] L.N. Trefethen, "Householder triangularization of a quasimatrix", IMA J
%   Numer Anal (2010) 30 (4): 887-897.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Developer note: Typically 'householder' will be slower than 'built-in', but
% should be more stable.
%
% Developer note: In the built-in approach, one could use a Chebyshev grid of
% double the size (rather than the Legendre grid), but given the complexity
% costs of QR vs. the cost of the transformation, it seems Legendre is
% preferable.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Deal with empty case:
if ( isempty(f) )
    R = [];
    E = [];
    return
end

% Default options:
% defaultMethod = 'built-in';
defaultMethod = 'householder';
defaultOutput = 'matrix';

if ( nargin < 3 || isempty(methodFlag) )
    methodFlag = defaultMethod;
end
if ( nargin < 2 || isempty(outputFlag) )
    outputFlag = defaultOutput;
end

% If f has only one column we simply scale it.
if ( size(f, 2) == 1 )
    R = sqrt(innerProduct(f, f));
    f = f./R;
    E = 1;
    return
end

% Simplify so that we don't do any extra work: (QR is O(m*n^2)? :/ )
f = simplify(f);

% We must enforce that f.values has at least as many rows as columns:
[n, m] = size(f);
if ( n < m )
    f = prolong(f, m);
    n = m;
end

% Decide which algorithm to use:
if ( strcmpi(methodFlag, 'householder') )
    % Call Trefethen's Householder implementation:
    [f, R, E] = qr_householder(f, outputFlag);
    return
else
    % The 'built-in' algorithm. I.e., qeighted discrete QR():
    if ( nargout == 3 )
        [f, R, E] = qr_builtin(f, outputFlag);
    else
        [f, R] = qr_builtin(f, outputFlag);
    end
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f, R, E] = qr_builtin(f, outputFlag)

% Project the values onto a Legendre grid: (where integrals of polynomials
% p_n*q_n will be computed exactly and on an n-point grid)
[n, m] = size(f);
xc = f.chebpts(n);
vc = f.barywts(n);
[xl, wl, vl] = legpts(n);
P = barymat(xl, xc, vc);
W = spdiags(sqrt(wl.'), 0, n, n);

% Compute the weighted QR factorisation:
if ( nargout == 3 )
    [Q, R, E] = qr(W * P * f.values, 0);
    % For consistency with the MATLAB QR behavior:
    if ( (nargin == 1) || ~(strcmpi(outputFlag, 'vector') || (outputFlag == 0)) )
        % Return E in matrix form:
        I = eye(m);
        E = I(:,E);
    end
else
    [Q, R] = qr(W * P * f.values, 0);
end

% Revert to the Chebyshev grid (and remove the weight and enforce diag(R) >= 0).
Winv = diag(1./sqrt(wl));   % Undo the weighting used for QR.
Pinv = barymat(xc, xl, vl); % Revert to Chebyshev grid (from Legendre).

% Enforce diag(R) >= 0.
s = sign(diag(R));
s(~s) = 1;
S = spdiags(s, 0, m, m);
Q = Pinv*Winv*Q*S;          % Fix Q.
R = S*R;                    % Fix R.

% Apply data to chebtech:
f.values = Q;                           % Adjust values of f.
f.coeffs = f.vals2coeffs(Q);            % Compute new coefficients.
f.vscale = max(abs(Q), [], 1);

% [TODO]: Update epslevel?

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f, R, E] = qr_householder(f, flag)
% See L.N. Trefethen, "Householder triangularization of a quasimatrix", IMA J
% Numer Anal (2010) 30 (4): 887-897. 

% Specify a tolerance:
tol = f.epslevel*max(f.vscale);

% Grab the size:
[n, m] = size(f);

% Create the Chebyshev nodes and quadrature weights of double the length:
% (Note: we double the size so that sum(f*f) is exact for f in P_{n}.)
x = f.chebpts(2*n);
w = f.quadwts(2*n);

    % Define the inner product as a nested function:
    function res = innerProd(a, b)
        res = w*(conj(a).*b);
    end

% Make the matrix of function values on a doubled grid so that Q*R=A:
A = get(prolong(f, 2*n), 'values');

% Generate the Legendre-Vandermonde matrix E via the standard recurrence:
E = ones(2*n, m);
E(:,2) = x;
for k = 3:m
    E(:,k) = ((2*k-3)*x.*E(:,k-1) - (k - 2)*E(:,k-2)) / (k - 1);
end
for k = 1:m
    E(:,k) = E(:,k) * sqrt(k - .5);
end

% Pre-allocate memory for R and V:
R = zeros(m);
V = zeros(2*n, m);

% Discretised version of code from Trefethen's paper:
for k = 1:m
    
    % Indices of the previous and following columns:
    I = 1:k-1; 
    J = k+1:m;
    scl = max(max(abs(E(:,k))), max(abs(A(:,k))));
    
    % Multiply the kth column of A with the basis in E:
    ex = innerProd(E(:,k), A(:,k));
    aex = abs(ex);
    
    % Adjust the sign of the kth column in E:
    if ( aex < eps*scl )
        s = 1; 
    else
        s = -sign(ex/aex);
    end
    E(:,k) = E(:,k) * s;
    
    % Compute the norm of the kth column of A:
    r = sqrt(innerProd(A(:,k), A(:,k)));
    R(k,k) = r;
    
    % Compute the reflection v:
    v = r*E(:,k) - A(:,k);
    % Make it more orthogonal:
    for j = I
        ev = innerProd(E(:,j), v);
        v = v - E(:,j)*ev;
    end
    
    % Normalize and store v:
    nv = sqrt(innerProd(v, v));
    if ( nv < tol*scl )
        v = E(:,k);
    else
        v = v / nv;
    end
    V(:,k) = v;
    
    % Subtract v from the remaining columns of A:
    for j = J
        av = innerProd(v, A(:,j));
        A(:,j) = A(:,j) - 2*v*av;
        rr = innerProd(E(:,k), A(:,j));
        A(:,j) = A(:,j) - E(:,k)*rr;
        R(k,j) = rr;
    end
    
end

% Form a discrete Q from the columns of V:
Q = E;
for k = m:-1:1
    for j = k:m
        vq = innerProd(V(:,k), Q(:,j));
        Q(:,j) = Q(:,j) - 2*V(:,k)*vq;
    end
end

% Compute the corresponding Chebyshev coefficients:
f.coeffs = f.vals2coeffs(Q);
% Trim the unneeded ones:
f.coeffs(1:n,:) = [];
% Comute new values:
f.values = f.coeffs2vals(f.coeffs);

% Update the vscale:
f.vscale = max(abs(Q), [], 1);

% Additional output argument:
if ( nargout == 3 )
    if ( nargin == 2 && strcmp(flag, 'vector') )
        E = 1:m;
    else
        E = eye(m);
    end
end

end
