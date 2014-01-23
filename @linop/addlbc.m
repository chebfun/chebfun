function L = addlbc(L, op, varargin)
%ADDLBC    Append a boundary condition at the left endpoint.
%   L = ADDLBC(L,VAL) sets a constraint that the function at the left
%   endpoint has value VAL.
%
%   L = ADDLBC(L,OP,VAL) imposes (OP*u) at the left endpoint equals VAL.
%
%   See also LINOP.ADDBC.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

if ( size(L,1) > 1 )
    error('ADDLBC syntax only applies to single-variable operators.')
end

% The domain of the LINOP
dom = L.domain;
% Create an feval FUNCTIONALBLOCK which allows us to impose the condition via
% addbc.
E = functionalBlock.feval(dom(1), dom);
if ( isnumeric(op) )
    % It's really just a boundary value.
    L = addbc(L, E, op);
else
    % Compose with the given operator. There is no special treatment given to
    % left/right boundary conditions, they're stored in the same way as general
    % boundary conditions imposed by calling addbc.
    L = addbc(L, E*op, varargin{:});
end

end
