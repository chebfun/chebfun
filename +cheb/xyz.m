function varargout = xyz
%XYZ   Three chebfun3 obejcts of the identity on [-1, 1, -1, 1, -1, 1].
%   CHEB.XYZ is shorthand for the expressions 
%   X = CHEBFUN3(@(X,Y,Z) X),
%   Y = CHEBFUN3(@(X,Y,Z) Y), and
%   Z = CHEBFUN3(@(X,Y,Z) Z).

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

x = chebfun3(@(x,y,z) x);
y = chebfun3(@(x,y,z) y);
z = chebfun3(@(x,y,z) z);

if nargout>1
    % For the syntax [x,y,z] = cheb.xyz
    varargout{1} = x;
    varargout{2} = y;
    varargout{3} = z;
else
    % For the syntax cheb.xyz we still want to put these variables into the
    % workspace:
    assignin('base', 'x', x);
    assignin('base', 'y', y);
    assignin('base', 'z', z);
end

end