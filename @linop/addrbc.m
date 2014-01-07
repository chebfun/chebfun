function L = addrbc(L, op, varargin)
%ADDRBC  Append to linop constraints (right BC)

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

if ( size(L,1) > 1 )
    error('ADDRBC syntax only applies to single-variable operators.')
end

d = L.domain;
E = linop.feval(d(end), d);
if isnumeric(op)
    % It's really just a boundary value.
    L = addbc(L, E, op);
else
    % Compose with the given operator.
    L = addbc(L, E*op, varargin{:});
end

end
