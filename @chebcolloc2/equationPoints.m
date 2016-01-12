function varargout = equationPoints(disc)
%EQUATIONPOINTS   Points at which collocation is enforced.
%   In CHEBCOLLOC2, functions are discretized at 2nd kind points but equations 
%   are enforced at 1st kind points, to avoid duplication at boundaries.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

pointsFun = @(n) chebtech1.chebpts(n);
[varargout{1:nargout}] = valsDiscretization.points(disc, pointsFun);

end
