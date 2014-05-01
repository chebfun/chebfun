function [Q, R, E] = qr(f, outputFlag, methodFlag)
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
%   to use in computing the QR factorisation. METHOD = 'built-in' will form a
%   weighted Legendre-Vandermonde matrix and orthogonalise this with the
%   standard Matlab QR algorithm. METHOD = 'householder' uses the technique
%   described in [1]. METHOD = 'householder' is the default option.
%
%   [1] L.N. Trefethen, "Householder triangularization of a quasimatrix", IMA J
%   Numer Anal (2010) 30 (4): 887-897.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
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
    Q = f;
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
    Q = f./R;
    E = 1;
    return
end

% TODO: This should probably be put back in if possible? NH Apr 2014
% Simplify so that we don't do any extra work: (QR is O(m*n^2)? :/ )
% f = simplify(f);

% Decide which algorithm to use:
if ( strcmpi(methodFlag, 'householder') )
    % Call Trefethen's Householder implementation:
    [Q, R, E] = qr_householder(f, outputFlag);
else
    % The 'built-in' algorithm. i.e., qeighted discrete QR():
    if ( nargout == 3 )
        [Q, R, E] = qr_builtin(f, outputFlag);
    else
        [Q, R] = qr_builtin(f, outputFlag);
    end
end

%% Update epslevel
% Since we don't know how to do this properly, we essentially assume that QR has
% condition number one. Therefore we assume Q has the same global accuracy as f,
% and simply factor out the new vscale. TODO: It may be sensible to include some
% knowledge of R here?
col_acc = f.epslevel.*f.vscale;  % Accuracy of each column in f.
glob_acc = max(col_acc);         % The best of these.
Q.epslevel = glob_acc./Q.vscale; % Scale out vscale of Q.

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f, R, E] = qr_builtin(f, outputFlag)

% We must enforce that f.coeffs has at least as many rows as columns:
[n, m] = size(f);
if ( n < m )
    f = prolong(f, m);
    n = m;
end

% Project the values onto a Legendre grid: (where integrals of polynomials
% p_n*q_n will be computed exactly and on an n-point grid)
[n, m] = size(f);
xc = f.chebpts(n);
vc = f.barywts(n);
[xl, wl, vl] = legpts(n);
P = barymat(xl, xc, vc);
W = spdiags(sqrt(wl.'), 0, n, n);

% Compute the weighted QR factorisation:
values = f.coeffs2vals(f.coeffs);
if ( nargout == 3 )
    [Q, R, E] = qr(W * P * values, 0);
    % For consistency with the MATLAB QR behavior:
    if ( (nargin == 1) || ...
        ~(strcmpi(outputFlag, 'vector') || isequal(outputFlag, 0)) )
        % Return E in matrix form:
        I = eye(m);
        E = I(:,E);
    end
else
    [Q, R] = qr(W * P * values, 0);
end

% Revert to the Chebyshev grid (and remove the weight and enforce diag(R) >= 0).
Winv = diag(1./sqrt(wl));   % Undo the weighting used for QR.
Pinv = barymat(xc, xl, vl); % Revert to Chebyshev grid (from Legendre grid).

% Enforce diag(R) >= 0.
s = sign(diag(R));
s(~s) = 1;
S = spdiags(s, 0, m, m);
Q = Pinv*Winv*Q*S;          % Fix Q.
R = S*R;                    % Fix R.

% Apply data to chebtech:
f.coeffs = f.vals2coeffs(Q);            % Compute new coefficients.
f.vscale = max(abs(Q), [], 1);

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f, R, Eperm] = qr_householder(f, flag)

% Get some useful values
[n, numCols] = size(f);
tol = max(f.epslevel.*f.vscale);

% Make the discrete analog of f:
newN = 2*max(n, numCols);
A = get(prolong(f, newN), 'values');

% Create the Chebyshev nodes and quadrature weights:
x = f.chebpts(newN);
w = f.quadwts(newN);

% Define the inner product as an anonymous function:
ip = @(f, g) w * (conj(f) .* g);

% Generate a discrete E (Legendre-Chebyshev-Vandermonde matrix) directly:
E = ones(size(A));
E(:,2) = x;
for k = 3:numCols % Recurrence relation:
    E(:,k) = ((2*k - 3)*x.*E(:,k-1) - (k - 2)*E(:,k-2)) / (k - 1);
end
% Scaling:
for k = 1:numCols
    E(:,k) = E(:,k) * sqrt((2*k - 1) / 2);
end
% Call the abstract QR method:
[Q, R] = abstractQR(A, E, ip, @(v) norm(v, inf), tol);

% Compute the corresponding Chebyshev coefficients:
f.coeffs = f.vals2coeffs(Q);
% Trim the unneeded ones:
f.coeffs(1:newN/2,:) = [];

% Update the vscale:
f.vscale = getvscl(f);

% Additional output argument:
if ( nargout == 3 )
    if ( nargin == 2 && strcmp(flag, 'vector') )
        Eperm = 1:numCols;
    else
        Eperm = eye(numCols);
    end
end

end
