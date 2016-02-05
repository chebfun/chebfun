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
%   described in [1]. METHOD = 'built-in' is the default option.
%
%   [1] L.N. Trefethen, "Householder triangularization of a quasimatrix", IMA J
%   Numer Anal (2010) 30 (4): 887-897.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
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
defaultMethod = 'built-in';
%defaultMethod = 'householder';
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

% Simplify so that we don't do any extra work: (QR is O(m*n^2)? :/ )
f = simplify(f);

% Decide which algorithm to use:
if ( strcmpi(methodFlag, 'householder') )
    % Call Trefethen's Householder implementation:
    [Q, R] = qr_householder(f);
else
    % The 'built-in' algorithm. i.e., qeighted discrete QR():
    [Q, R] = qr_builtin(f);
end

% Additional output argument:
if ( nargout == 3 )
    [~, numCols] = size(f);
    if ( nargin > 1 && strcmp(outputFlag, 'vector') )
        E = 1:numCols;
    else
        E = eye(numCols);
    end
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f, R] = qr_builtin(f)

% We must enforce that f.coeffs has at least as many rows as columns:
[n, m] = size(f);
if ( n < m )
    f = prolong(f, m);
end

n = size(f.coeffs, 1); 
scl = sqrt((0:n-1)'+.5);
legA = cheb2leg( f.coeffs );          % Get Legendre coefficients.
legA = bsxfun(@times, legA, 1./scl ); % Orthonormalization constants.
[Qleg, R] = qr(legA, 0);              % Discrete QR
     
% Ensure diagonals of R are of +ve sign.
s = sign(diag(R)); 
s(~s) = 1; 
S = spdiags(s, 0, m, m); 
R = S*R; 
Qleg = Qleg*S;  

Qleg = bsxfun(@times, Qleg, scl);    % Apply scaling.
Qcheb = leg2cheb( Qleg );            % Convert back to Chebyshev.
f.coeffs = Qcheb;                    % Store coefficients. 

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f, R] = qr_householder(f)

% Get some useful values
[n, numCols] = size(f);
tol = eps*max(vscale(f));

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
f.coeffs(newN/2+1:end,:) = [];

end
