function L = addrbc(L, op, varargin)
%ADDLBC    Append a boundary condition at the right endpoint.
%   L = ADDRBC(L,VAL) sets a constraint that the function at the right
%   endpoint has value VAL.
%
%   L = ADDLBC(L,OP,VAL) imposes (OP*u) at the right endpoint equals VAL.
%
%   See also LINOP.ADDBC.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

if ( size(L,1) > 1 )
    error('ADDRBC syntax only applies to single-variable operators.')
end

d = L.domain;
E = functionalBlock.feval(d(end), d);
if isnumeric(op)
    % It's really just a boundary value.
    L = addbc(L, E, op);
else
    % Compose with the given operator.
    L = addbc(L, E*op, varargin{:});
end

end
