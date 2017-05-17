function f = compose(f, op)
%COMPOSE   Compose command for DISKFUNV objects.
%   F = COMPOSE(F, G) returns the composition G(F) of the DISKFUNV object F
%   and a CHEBFUN2 or CHEBFUN2V with two components G.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty CHEBFUN objects:
if ( isempty(f) || isempty(op) )
    if ( isa(op, 'chebfun2v') )
        f = diskfunv();
    else
        f = diskfun();
    end
end

% F is real and has two components, so we can compose with a CHEBFUN2 or
% CHEBFUN2V.

if ( isa(op, 'chebfun2') )
    % Get components of F:
    f1 = f.components{1};
    f2 = f.components{2};
    
    % Check that image(f) is in domain(op):
    rangef = minandmax2est(f);              % Estimate of image(f).
    tol = 100 * chebfun2eps * max(vscale(f), vscale(op)) * ...
            norm(f1.domain, inf);           % Tolerance.
    if ( ~isSubset(rangef, op.domain, tol) )
        error('CHEBFUN:DISKFUNV:COMPOSE:DomainMismatch2', ...
            'OP(F) is not defined, since image(F) is not contained in domain(OP).')
    end
    
    % Call constructor:
    f = diskfun(@(x,y) op(feval(f1, x, y), feval(f2, x, y)));
    
elseif ( isa(op, 'chebfun2v') )
    % Check that OP has two components:
    if ( op.nComponents ~= 2 )
        error('CHEBFUN:DISKFUNV:compose:nComponents', ...
            'CHEBFUN2V(DISKFUNV) is defined for CHEBFUN2V objects with 2 components.')
    end
    
    % Get components of op:
    op1 = op(1);
    op2 = op(2);
    
    % Get components of f:
    f1 = f.components{1};
    f2 = f.components{2};
    
    % Check that image(f) is in domain(op):
    rangef = minandmax2est(f);              % Estimate of image(f).
    tol = 100 * chebfun2eps * max(vscale(f), vscale(op)) * ...
            norm(f1.domain, inf);           % Tolerance.
    if ( ~isSubset(rangef, op1.domain, tol) )
        error('CHEBFUN:DISKFUNV:COMPOSE:DomainMismatch2v', ...
            'OP(F) is not defined, since image(F) is not contained in domain(OP).')
    end
    
    % Call constructor:
    f = diskfunv(@(x,y) op1(feval(f1, x, y), feval(f2, x, y)), ...
        @(x,y) op2(feval(f1, x, y), feval(f2, x, y)));
    
else
    error('CHEBFUN:DISKFUNV:COMPOSE:OP', ...
        'Can compose a DISKFUNV with a CHEBFUN2 or CHEBFUN2V.')
end

end
