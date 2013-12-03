function L = addlbc(L, op, varargin)
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
%ADDLBC  Append to linop constraints (left BC)
d = L.domain;
E = linop.feval(d(1), d);
L = addbc(L, E*op, varargin{:});
end
