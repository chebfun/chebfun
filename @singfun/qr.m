function [f, R, E] = qr(f, outputFlag, methodFlag)
%QR   QR factorisation of a SINGFUN.
%   [Q, R] = QR(F) returns a QR factorisation of F such that F = Q*R, where the
%   SINGFUN Q is normalized (with respect to the continuous L^2 norm on [-1,1])
%   and has only a single column. R is a positive scalar which is square root of
%   the L^2 norm of F. Note that Q has only one column since SINGFUN does not
%   support array-valued SINGFUN objects.
%
%   [Q, R, E] = QR(F) produces the same Q and R, along with a trivial 
%   permutation matrix E, which is just the scalar 1.
%
%   [Q, R, E] = QR(F, 'vector') or [Q, R, E] = QR(F, 'matrix') returns exactly 
%   the same result as the second input argument is omitted.
%
%   QR(F, 'vector', METHOD) or QR(F, 'vector', METHOD) allows a third input 
%   which is expected to be a string. This is only to keep the calling sequence 
%   to be consistent with that of the SMOOTHFUN counterpart.
%
% See also innerProduct.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Simply call innerProduct:
R = sqrt(innerProduct(f, f));
if ( isinf(R) )
    error('CHEBFUN:SINGFUN:qr:infNorm', ...
        'The L^2 norm of the input SINGFUN is infinite.')
end

% Normalization:
f = f/R;

% The permutation scalar is alway 1:
E = 1;

end
