function L = addrbc(L, op, varargin)
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
%ADDRBC  Append to linop constraints (right BC)
d = L.domain;
E = linop.feval(d(end), d);
L = addbc(L, E*op, varargin{:});
end
