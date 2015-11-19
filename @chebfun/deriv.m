function out = deriv(f, xx, m)
%DERIV   Evaluate a derivative of a CHEBFUN.
%
%   OUT = DERIV(F, X) evaluates the first derivative of the CHEBFUN F at the
%   points in X. If F is a quasimatrix with columns F1, ..., FN, then the result
%   will be [F1(X), ..., FN(X)], the horizontal concatenation of the results of
%   evaluating each column at the points in X.
%
%   OUT = DERIV(F, X, M) is the same as above, but returns the values of the Mth
%   derivative of F.
%
%   Example:
%     u = chebfun(@sin);
%     deriv(u, 1) % u'(1)
%     deriv(u, .5, 2) % u''(.5)
%     deriv(u, 0.1:0.01:0.2) % u'(0.1:0.01:0.2)
%
% See also chebfun/diff, chebfun/feval.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% The most trivial case:
if ( isempty(f) )
    return
end

if ( nargin < 3 )
    m = 1; % By default, compute first derivative
end

if ( m == 0 )
    % Trivial case
    out = feval(f, xx);
else
    out = feval(diff(f, m), xx);
end

end
