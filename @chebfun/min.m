function [y, x] = min(f, varargin)
%MIN    Minimum values of a chebfun.
%   MIN(F) returns the minimum value of the chebfun F.
%
%   [Y, X] = MINANDMAX(F) returns also point X such that F(X) = Y.
%
%   [Y, X] = MINANDMAX(F, 'local') returns not just the global minimum value,
%   but all of the local minima.
%
%   If F is complex-valued, absolute values are taken to determine the minima,
%   but the resulting values correspond to those of the original function.
%
% See also CHEBFUN/MAX, CHEBFUN/MINANDMAX.
%
% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call MINANDMAX():
[y, x] = minandmax(f, varargin{:});

% Extract the minimum:
y = y(:,1);
x = x(:,1);

% If the 'local' flag is passed, we must do some work to determine which are
% minima:
if ( numel(y) > 1 )
    % For interior extrema, look at 2nd derivative:
    idx = feval(diff(f, 2), x) > 0;
    % For end-points, look at 1st derivative:
    idx([1, end]) = feval(diff(f), x([1, end])).*[1, -1]' > 0;
    % Extract relelvant entries:
    y = y(idx);
    x = x(idx);
end

end