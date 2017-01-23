function varargout = xyzsphere
%XYDISK  Diskfun objects for x and y on the disk.
%   [X, Y] = CHEB.XYDISK returns diskfun objects for the functions 
%   @(x,y) x, and  @(x,y) y defined on the unit disk.
%
%   CHEB.XYDISK is shorthand for the expressions 
%   X = DISKFUN(@(X,Y) X) and
%   Y = DISKFUN(@(X,Y) Y). 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

x = diskfun(@(x,y) x);
y = diskfun(@(x,y) y);

if ( nargout > 0 )
    
    % For the syntax [x, y, z] = cheb.xyzsphere:
    varargout{1} = x;
    varargout{2} = y;
    
    if ( nargout > 2 ) 
        error('CHEB:XYDISK:TooManyOutputs',... 
            'Too many output arguments. CHEB.XYDISK only returns "x", and "y".')
    end
    
else
    % Put 'x', 'y', and 'z' into the workspace:
    assignin('base', 'x', x)
    assignin('base', 'y', y)

end
end
