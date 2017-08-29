function f = compose(f, op)
%COMPOSE   Compose command for SPHEREFUNV objects.
%   F = COMPOSE(F, G) returns the composition G(F) of the SPHEREFUNV object
%   F and a CHEBFUN3 or CHEBFUN3V with three components G.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty CHEBFUN objects:
if ( isempty(f) || isempty(op) )
    if ( isa(op, 'chebfun3v') )
        f = spherefunv();
    else
        f = spherefun();
    end
end

% F is real and has three components, so we can compose with a CHEBFUN3 or
% CHEBFUN3V.
% TODO? compose when op is a spherefun object?

if ( isa(op, 'chebfun3') )
    % Get components of F:
    f1 = f.components{1};
    f2 = f.components{2};
    f3 = f.components{3};
    
    % Check that image(f) is in domain(op):
    rangef = minandmax2est(f);              % Estimate of image(f).
    tol = 100 * chebfun2eps * max(vscale(f), vscale(op)) * ...
            norm(f1.domain, inf);           % Tolerance.
    if ( ~isSubset(rangef, op.domain, tol) )
        error('CHEBFUN:SPHEREFUNV:COMPOSE:DomainMismatch3', ...
            'OP(F) is not defined, since image(F) is not contained in domain(OP).')
    end
    
    % Call constructor:
    f = spherefun(@(x,y,z) op(feval(f1, x, y, z), feval(f2, x, y, z), ...
        feval(f3, x, y, z)));
    
elseif ( isa(op, 'chebfun3v') )
    % Check that OP has three components:
    if ( op.nComponents ~= 3 )
        error('CHEBFUN:SPHEREFUNV:compose:nComponents', ...
            'CHEBFUN3V(SPHEREFUNV) is defined for CHEBFUN3V objects with 3 components.')
    end
    
    % Get components of op:
    op1 = op(1);
    op2 = op(2);
    op3 = op(3);
    
    % Get components of f:
    f1 = f.components{1};
    f2 = f.components{2};
    f3 = f.components{3};
    
    % Check that image(f) is in domain(op):
    rangef = minandmax2est(f);              % Estimate of image(f).
    tol = 100 * chebfun2eps * max(vscale(f), vscale(op)) * ...
            norm(f1.domain, inf);           % Tolerance.
    if ( ~isSubset(rangef, op1.domain, tol) )
        error('CHEBFUN:SPHEREFUNV:COMPOSE:DomainMismatch3v', ...
            'OP(F) is not defined, since image(F) is not contained in domain(OP).')
    end
    
    % Calling directly the constructor for [ op1(f), op2(f), op3(f) ] as below
    % seems to be slightly faster than building a spherefun for each component
    % and then calling the spherefunv constructor on the three spherefuns.
    
    % Call constructor:
    f = spherefunv(@(x,y,z) op1(feval(f1, x, y, z), feval(f2, x, y, z), ...
        feval(f3, x, y, z)), ...
        @(x,y,z) op2(feval(f1, x, y, z), feval(f2, x, y, z), ...
        feval(f3, x, y, z)), ...
        @(x,y,z) op3(feval(f1, x, y, z), feval(f2, x, y, z), ...
        feval(f3, x, y, z)));
    
else
    error('CHEBFUN:SPHEREFUNV:COMPOSE:OP', ...
        'Can compose a SPHEREFUNV with a CHEBFUN3 or CHEBFUN3V.')
end

end