function [Q, R, E] = qr(f, outputFlag, methodFlag)
%QR   QR factorisation of an array-valued TRIGTECH.
%   [Q, R] = QR(F) returns a QR factorisation of F such that F = Q*R, where the
%   TRIGTECH Q is orthogonal (with respect to the continuous L^2 norm on [-1,1])
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
%   QR(F, 'vector', METHOD) or QR(F, 'vector', METHOD) specifies which
%   method to use in computing the QR factorisation. METHOD = 'built-in'
%   will form a weighted Vandermonde matrix and orthogonalise this with the
%   standard Matlab QR algorithm. METHOD = 'householder' uses the technique
%   described in [1]. METHOD = 'built-in' is the default option.
%
%   [1] L.N. Trefethen, "Householder triangularization of a quasimatrix", IMA J
%   Numer Anal (2010) 30 (4): 887-897.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case:
if ( isempty(f) )
    Q = f;
    R = [];
    E = [];
    return
end

% Default options:
defaultMethod = 'built-in';
% defaultMethod = 'householder';
defaultOutput = 'matrix';

if ( nargin < 3 || isempty(methodFlag) )
    methodFlag = defaultMethod;
end
if ( nargin < 2 || isempty(outputFlag) )
    outputFlag = defaultOutput;
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

end

function [f, R, Eperm] = qr_householder(f, flag)

% Get some useful values
[n, numCols] = size(f);
tol = max(eps*vscale(f));

% Make the discrete analog of f:
newN = 2*max(n, numCols);
A = get(prolong(f, newN), 'values');

% Create the trigonometric nodes and quadrature weights:
x = f.trigpts(newN);
w = f.quadwts(newN);

% Define the inner product as an anonymous function:
ip = @(f, g) w * (conj(f) .* g);

% Work with sines and cosines instead of complex exponentials.
E1 = cos(pi*x*(0:floor(numCols/2))); E1(:,1) = E1(:,1)/sqrt(2); 
E2 = sin(pi*x*(1:ceil(numCols/2)-1));
E = zeros(size(A));
E(:,[1 2:2:end]) = E1; 
E(:,3:2:end) = E2; 

% Call the abstract QR method:
[Q, R] = abstractQR(A, E, ip, @(v) norm(v, inf), tol);

f.values = Q;
f.coeffs = f.vals2coeffs(Q); 
% If any columns of f were not real, we cannot guarantee that the columns
% of Q should remain real.
f.isReal(:) = all(f.isReal);

% Prune the unneeded coefficients.
f = prolong(f,n);

% Additional output argument:
if ( nargout == 3 )
    if ( nargin == 2 && strcmp(flag, 'vector') )
        Eperm = 1:numCols;
    else
        Eperm = eye(numCols);
    end
end

end

function [f, R, E] = qr_builtin(f, outputFlag)
% Quadratures are being done using trapezoidal rule. We must enforce that 
% f.coeffs has at least as many rows as columns
[nf, mf] = size(f);
n = max(nf, mf);        % We want Q to have mf columns. 
isReal = all(f.isReal);
f = prolong(f, n);
if (isReal)
    % Force values to be real if f is real (avoid complex rounding errors in R)
    f.values = real(f.values);
end

if ( nargout == 3 )
    [Q, R, E] = qr(f.values, 0);
    % For consistency with the MATLAB QR behavior:
    if ( ~(strcmpi(outputFlag, 'vector') || isequal(outputFlag, 0)) )
        % Return E in matrix form:
        I = eye(mf);
        E = I(:,E);
    end
else
    [Q, R] = qr(f.values, 0);
end

% Enforce diag(R) >= 0:
s = sign(diag(R));
s(~s) = 1;
S = spdiags(s, 0, mf, mf);
Q = Q*S; % Fix Q.
R = S*R; % Fix R.

% Scaling:
w = 2/n; % Weights in trapezoidal rule.
W = sqrt(w);
Q = Q/W; % Scale Q.
R = W*R; % Scale R.

% Apply data to TRIGTECH:
f.values = Q;
f.coeffs = f.vals2coeffs(Q); 
% If any columns of f were not real, we cannot guarantee that the columns
% of Q should remain real.
f.isReal(:) = isReal;
% Prune the unneeded coefficients.
f = prolong(f, nf);

end

