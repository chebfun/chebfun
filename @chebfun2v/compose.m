function f = compose(f, op)
%COMPOSE   Compose command for CHEBFUN2V objects.
%   F = COMPOSE(F, G) returns the composition G(F) of the real-valued CHEBFUN2V
%   object F and a suitable CHEBFUN* object G (* may be 2, 2V, 3, or 3V).

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TO DO: Add test that range(f) is in domain(OP).

% Deal with empty CHEBFUN objects:
if ( isempty(f) || isempty(op) )
    if ( isa(op, 'chebfun2v') || isa(op, 'chebfun3v') )
        f = chebfun2v();
    else
        f = chebfun2();
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
    error('CHEBFUN:CHEBFUN2V:COMPOSE:Complex', ...
        'The first CHEBFUN2V object must be real-valued.')
end


if ( f.nComponents == 1 )
    % If F has only one component, send to CHEBFUN2 compose:
    f = compose(f.components{1}, op);
    
elseif ( f.nComponents == 2 )
    % f has 2 components, so we can compose with a CHEBFUN2 or CHEBFUN2V.
    
    if ( isa(op, 'chebfun2') )
        f1 = f.components{1};
        f2 = f.components{2};
        
        % Call constructor:
        f = chebfun2(@(x,y) op(feval(f1, x, y), feval(f2, x, y)), f1.domain);
        
    elseif ( isa(op, 'chebfun2v') )
        % Call compose for each component:
        F = compose(f, op.components{1});
        for jj = 2:op.nComponents
            F = [F; compose(f, op.components{jj})];
        end
        f = F;
        
    else
        error('CHEBFUN:CHEBFUN2V:COMPOSE:OP2', 'OP not recognized.')
        
    end
    
elseif ( f.nComponents == 3 )
    % f has 3 components, so we can compose with a CHEBFUN3 or CHEBFUN3V.
    
    if ( isa(op, 'chebfun3') )
        % Get components:
        f1 = f.components{1};
        f2 = f.components{2};
        f3 = f.components{3};
        
        % Call constructor:
        f = chebfun2(@(x,y) op(feval(f1, x, y), feval(f2, x, y), ...
            feval(f3, x, y)), f1.domain);
        
    elseif ( isa(op, 'chebfun3v') )
        % Call compose for each component:
        F = compose(f, op(1));
        for jj = 2:op.nComponents
            F = [F; compose(f, op(jj))];
        end
        f = F;
        
    else
        error('CHEBFUN:CHEBFUN2V:COMPOSE:OP3', 'OP not recognized.')
        
    end
    
end

end
