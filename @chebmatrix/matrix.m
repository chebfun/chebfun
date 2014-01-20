function M = matrix(A,varargin)

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

d = A.prefs.discretization(A,varargin{:});
M = matrix(d);

end