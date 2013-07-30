function [y, x] = max(f, varargin)
%MAX    Maximum values of a chebfun.
%   MAX(F) returns the maximum value of the chebfun F.
%
%   [Y, X] = MAX(F) returns also point X such that F(X) = Y.
%
%   [Y, X] = MAX(F, 'local') returns not just the global maximum value,
%   but all of the local maxima.
%
%   If F is complex-valued, absolute values are taken to determine the maxima,
%   but the resulting values correspond to those of the original function.
%
% See also CHEBFUN/MIN, CHEBFUN/MINANDMAX.
%
% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call MINANDMAX():
[y, x] = minandmax(f, varargin{:});

% Extract the maximum:
y = y(:,2);
x = x(:,2);

% If the 'local' flag is passed, we must do some work to determine which are
% maxima:
if ( numel(y) > 1 )
    % For interior extrema, look at 2nd derivative:
    idx = feval(diff(f, 2), x) > 0;
    % For end-points, look at 1st derivative:
    idx([1, end]) = feval(diff(f), x([1, end])).*[1, -1]' < 0;
    % Extract relelvant entries:
    y = y(idx);
    x = x(idx);
end

end