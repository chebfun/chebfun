function [Q, R, E] = qr(f, varargin)
%QR   QR factorisation of an array-valued BNDFUN.
%   [Q, R] = QR(F) returns a QR factorisation of F such that F = Q*R, where the
%   BNDFUN Q is orthogonal (with respect to the continuous L^2 norm on the
%   domain of F) and of the same size as F and R is an m x m upper-triangular
%   matrix when F has m columns.
%
%   [Q, R, E] = QR(F) produces unitary Q, upper-triangular R, and a permutation
%   matrix E so that F*E = Q*R. The column permutation E is chosen to reduce
%   fill-in in R.
%
%   [Q, R, E] = QR(F, 'vector') returns the permutation information as a vector
%   instead of a matrix.  That is, E is a row vector such that F(:,E) = Q*R.
%   Similarly, [Q, R, E] = QR(F, 'matrix') returns a permutation matrix E. This
%   is the default behavior.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case:
if ( isempty(f) )
    Q = [];
    R = [];
    E = [];
    return
end

% Initialise Q to be a BNDFUN:
Q = f;

% Rescaling factor, (b - a)/2.
rescaleFactor = .5*diff(f.domain);

% Call QR on the ONEFUN of f:
if ( nargout == 3 )
    [Q.onefun, R, E] = qr(f.onefun, varargin{:});
else
    [Q.onefun, R] = qr(f.onefun, varargin{:});
end

% Rescale so that columns of Q will be orthonormal (rather than orthogonal):
Q = Q/sqrt(rescaleFactor);

% Rescale R so that f = QR:
R = R*sqrt(rescaleFactor);

end
