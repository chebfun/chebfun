function v = grad(f)
%GRAD Gradient of a BALLFUN in cartesian coordinates.
%   GRAD(F) is the gradient of the BALLFUN F expressed in
%   cartesian coordinates.
%
%   This is shorthand for the command GRADIENT.
%
% See also DIV, CURL

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

v = gradient(f);
end