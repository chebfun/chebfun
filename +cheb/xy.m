function varargout = xy
%XY   Chebfun2 objects for the @(x,y)x and @(x,y)y on [-1,1]^2.
%   [X, Y] = CHEB.XY returns chebfun2 objects for the functions @(x,y)x 
%   and @(x,y)y defined on [-1,1]^2.
%   
%   CHEB.XY is shorthand for the expressions 
%   CHEBFUN2(@(X,Y) X), and 
%   CHEBFUN2(@(X,Y) Y).

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

x = chebfun2(@(x,y) x);
y = chebfun2(@(x,y) y);

if ( nargout > 0 )
    
    % For the syntax [x, y] = cheb.xy:
    varargout{1} = x;
    varargout{2} = y;
    
    if ( nargout > 2 ) 
        error('CHEB:XY:TooManyOutputs',... 
            'Too many output arguments. CHEB.XY only returns "x" and "y".')
    end
    
else
    % Put 'x' and 'y' into the workspace:
    assignin('base', 'x', x)
    assignin('base', 'y', y)
    
    % And display the variables:
    display(x)
    display(y)
end
end