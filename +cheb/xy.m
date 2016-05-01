function varargout = xy
%XY   A chebfun2 of the @(x,y)x and @(x,y)y on [-1,1]^2.
%   CHEB.XY is shorthand for the expressions 
%   CHEBFUN2(@(X,Y) X), and 
%   CHEBFUN2(@(X,Y) Y).

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

x = chebfun2(@(x,y) x);
y = chebfun2(@(x,y) y);

if ( nargout > 1 )
    % For the syntax [x, y] = cheb.xy:
    varargout{1} = x;
    varargout{2} = y;
else
    % For the syntax cheb.xy we still want to put these variables into the 
    % workspace:
    assignin('base', 'x', x)
    assignin('base', 'y', y)
    % and display
    display(x)
    display(y)
end

end