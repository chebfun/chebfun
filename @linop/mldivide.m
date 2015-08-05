function u = mldivide(L, f, varargin)
%\         Solve a linear system (same as LINOP.LINSOLVE).
%   Important: A CHEBOPPREF object PREFS has to be passed. When this method
%   is called via CHEBOP/MLDIVIDE, PREFS is inherited from the CHEBOP level.
%
%  Copyright 2015 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.
u = linsolve(L, f, varargin{:});

end
