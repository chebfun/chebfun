function varargout = equationPoints(disc)
%EQUATIONPOINTS   Points at which TRIGCOLLOC discretization is enforced.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

pointsFun = @(n) trigtech.trigpts(n);
[varargout{1:nargout}] = valsDiscretization.points(disc, pointsFun);

end
