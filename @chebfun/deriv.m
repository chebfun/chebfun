function out = deriv(f, xx, varargin)
%DERIV   Evaluate a derivative of a CHEBFUN.
%
%   DERIV(F, X) evaluates the first derivative of the CHEBFUN F at the points in
%   X. If F is a quasimatrix with columns F1, ..., FN, then the result will be
%   [F1(X), ..., FN(X)], the horizontal concatenation of the results of
%   evaluating each column at the points in X.
%
%   DERIV(F, X, M) is the same as above, but returns the values of the Mth
%   derivative of F.
%
%   DERIV(F, X, S) or DERIV(F, X, S, M) where S is one of the strings 'left',
%   'right', '+', or '-', evaluates the left- or right-sided limit as described
%   in chebfun/feval.
%
%   Example 1:
%     u = chebfun(@sin);
%     deriv(u, 1) % u'(1)
%     deriv(u, .5, 2) % u''(.5)
%     deriv(u, 0.1:0.01:0.2) % u'(0.1:0.01:0.2)
%
%   Example 2:
%     x = chebfun('x')
%     u = abs(x);
%     deriv(u, 0, 'left')  % ans = -1
%     deriv(u, 0, 'right') % ans = +1
%
% See also chebfun/diff, chebfun/feval.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% The most trivial case:
if ( isempty(f) )
    return
end

% By default, compute first derivative:
m = 1;      

% Parse the inputs:
if ( nargin == 3 )
    if ( isnumeric(varargin{1}) ) % DERIV(F, X, M)
        m = varargin{1};
        varargin = {};
    else                          % DERIV(F, X, S)
        m = 1;  % By default, compute first derivative
    end
elseif ( nargin == 4 )            % DERIV(F, X, S, M)
    m = varargin{2};
    varargin{2} = [];
end

if ( m == 0 )
    % Trivial case
    out = feval(f, xx, varargin{:});
else
    out = feval(diff(f, m), xx, varargin{:});
end

end
