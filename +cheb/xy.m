function [out1, out2] = xy
%XY   Two chebfun2 obejcts of the identity on [-1, 1, -1, 1].
%   CHEB.XY is shorthand for the expressions 
%   X = CHEBFUN2(@(X,Y) X) and
%   Y = CHEBFUN2(@(X,Y) Y).

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

out1 = chebfun2(@(x,y) x);
out2 = chebfun2(@(x,y) y);

end