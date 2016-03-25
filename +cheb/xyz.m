function [out1, out2, out3] = xyz
%XYZ   Three chebfun3 obejcts of the identity on [-1, 1, -1, 1, -1, 1].
%   CHEB.XYZ is shorthand for the expressions 
%   X = CHEBFUN3(@(X,Y,Z) X),
%   Y = CHEBFUN3(@(X,Y,Z) Y), and
%   Z = CHEBFUN3(@(X,Y,Z) Z).

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

out1 = chebfun3(@(x,y,z) x);
out2 = chebfun3(@(x,y,z) y);
out3 = chebfun3(@(x,y,z) z);

end