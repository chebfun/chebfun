function u = mldivide(L, f, varargin)
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
u = linsolve(L, f, varargin{:});
end
