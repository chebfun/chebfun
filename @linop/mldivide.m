function u = mldivide(L, f, varargin)
%\         Solve a linear system (same as LINOP.LINSOLVE).

%  Copyright 2014 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.
u = linsolve(L, f, varargin{:});

end
