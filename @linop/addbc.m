function L = addbc(L, varargin)
%ADDBC  Append to linop constraints (keep existing).

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

L.constraint = append(L.constraint, varargin{:});

end
