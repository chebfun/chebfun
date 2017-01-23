function f = compose(f, op)
%COMPOSE   Compose command for CHEBFUN3V objects.
%   F = COMPOSE(F, G) returns the composition G(F) of the real-valued CHEBFUN3V
%   object F and a suitable CHEBFUN* object G (* may be 2, 2V, 3, or 3V).

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty CHEBFUN objects:
if ( isempty(f) || isempty(op) )
    if ( isa(op, 'chebfun2v') || isa(op, 'chebfun3v') )
        f = chebfun3v();
    else
        f = chebfun3();
    end
end

% If OP is a CHEBFUN2 or CHEBFUN2v object, F may be scalar and complex-valued.
if ( ( f.nComponents == 1 ) && ( isa(op, 'chebfun2') || isa(op, 'chebfun2v') ) )
    % Interpret OP(f) as OP(real(f), imag(f))
    F = [real(f); imag(f)];
    f = compose(F, op);
    return
end

% In any other case f must be real.
if ( ~isreal(f) )
    error('CHEBFUN:CHEBFUN3V:COMPOSE:Complex', ...
        'The first CHEBFUN3V object must be real-valued.')
end


if ( f.nComponents == 1 )
    % If F has only one component, send to CHEBFUN3 compose:
    f = compose(f.components{1}, op);
    
elseif ( f.nComponents == 2 )
    % f has 2 components, so we can compose with a CHEBFUN2 or CHEBFUN2V.
    
    if ( isa(op, 'chebfun2') )
        % Get components:
        f1 = f.components{1};
        f2 = f.components{2};
        
        % Check that the image of f is in the domain of op:
        rangef = minandmax3est(f);          % Estimate of image(f).
        tol = 100 * chebfun3eps * max(vscale(f), vscale(op)) * ...
            norm(f1.domain, inf);           % Tolerance.
        if ( ~isSubset(rangef, op.domain, tol) )
            error('CHEBFUN:CHEBFUN3V:COMPOSE:DomainMismatch2', ...
                'OP(F) is not defined, since image(F) is not contained in domain(OP).')
        end
        
        % If f is periodic, then OP(f) should be periodic:
        pref = chebfunpref;
        if ( isPeriodicTech(f) )
            pref.tech = get(f1.rows(:,1).funs{1}, 'tech');
        end
        
        % Call constructor:
        f = chebfun3(@(x,y,z) op(feval(f1, x, y, z), feval(f2, x, y, z)), ...
            f1.domain, pref);
        
    elseif ( isa(op, 'chebfun2v') )
        F = compose(f, op(1));
        for jj = 2:op.nComponents
            F = [F; compose(f, op(jj))];
        end
        f = F;
        
    else
        error('CHEBFUN:CHEBFUN3V:COMPOSE:OP2', ...
            'Can compose only with a CHEBFUN2 or CHEBFUN2V, since F has 2 components.')
        
    end
    
elseif ( f.nComponents == 3 )
    % f has 3 components, so we can compose with a CHEBFUN3 or CHEBFUN3V.
    
    if ( isa(op, 'chebfun3') )
        % Get components:
        f1 = f.components{1};
        f2 = f.components{2};
        f3 = f.components{3};
        
        % Check that the image of f is in the domain of op:
        rangef = minandmax3est(f);          % Estimate of image(f).
        tol = 100 * chebfun3eps * max(vscale(f), vscale(op)) * ...
            norm(f1.domain, inf);           % Tolerance.
        if ( ~isSubset(rangef, op.domain, tol) )
            error('CHEBFUN:CHEBFUN3V:COMPOSE:DomainMismatch3', ...
                'OP(F) is not defined, since image(F) is not contained in domain(OP).')
        end
        
        % If f is periodic, then OP(f) should be periodic:
        pref = chebfunpref;
        if ( isPeriodicTech(f) )
            pref.tech = get(f1.rows(:,1).funs{1}, 'tech');
        end
        
        % Call constructor:
        f = chebfun3(@(x,y,z) op(feval(f1, x, y, z), feval(f2, x, y, z), ...
            feval(f3, x, y, z)), f1.domain, pref);
        
    elseif ( isa(op, 'chebfun3v') )
        F = compose(f, op.components{1});
        for jj = 2:op.nComponents
            F = [F; compose(f, op.components{jj})];
        end
        f = F;
        
    else
        error('CHEBFUN:CHEBFUN3V:COMPOSE:OP3', ...
            'Can compose only with a CHEBFUN3 or CHEBFUN3V, since F has 3 components.')
        
    end
    
end

end
