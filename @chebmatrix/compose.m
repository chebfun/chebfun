function f = compose(f, op)
%COMPOSE   Compose command for CHEBMATRIX objects.
%   H = COMPOSE(F, G) with a CHEBFUN2 G returns H(t) = G(F1(t), F2(t)), where F
%   is a 2 by 1 CHEBMATRIX whose entries F1, F2 are real-valued CHEBFUNs.
%   If G is a CHEBFUN2V, H is a CHEBMATRIX.
%
%   H = COMPOSE(F, G) with a CHEBFUN3 G returns H(t) = G(F1(t), F2(t), F3(t)),
%   where F is a 3 by 1 CHEBMATRIX whose entries F1, F2, F3 are real-valued
%   CHEBFUNs. If G is a CHEBFUN3V, H is a CHEBMATRIX.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TO DO: Add test that range(f) is in domain(OP).

% Deal with empty CHEBFUN objects:
if ( isempty(f) || isempty(op) )
    f = chebfun();
end

[m, n] = size(f);
% Ensure that n = 1.
if ( n ~= 1 )
    error('CHEBFUN:CHEBMATRX:compose:columns', ...
        'Can compose only with an 2 by 1 or 3 by 1 chebmatrix.')
end

% Only m = 2 and m = 3 are supported.
if ( m == 2 )
    % Can compose with a CHEBFUN2 or CHEBFUN2V.
    
    % Extract components:
    x = f.blocks{1};
    y = f.blocks{2};
    
    % Make sure all components are real-valued.
    if ( ~isreal(x) || ~isreal(y) )
        error('CHEBFUN:CHEBMATRIX:compose:complex2', ...
            'Can compose only with a CHEBMATRIX with real-valued CHEBFUN entries.')
    end
    
    if ( isa(op, 'chebfun2') )
        f = chebfun(@(t) op(feval(x, t), feval(y, t)), x.domain);
        
    elseif ( isa(op, 'chebfun2v') )
        % Call compose for each component of OP.
        F = compose(f, op.components{1});
        for jj = 2:op.nComponents
            F = [ F; compose(f, op(jj)) ];  % op.components{jj} does not work.
        end
        f = F;
        
    else
        error('CHEBFUN:CHEBMATRIX:compose:entries2', ...
            'Can compose a 2 by 1 CHEBMATRIX only with a CHEBFUN2 or CHEBFUN2V.')
    end
    
elseif ( m == 3 )
    % Can compose with a CHEBFUN3 or CHEBFUN3V
    
    % Extract components:
    x = f.blocks{1};
    y = f.blocks{2};
    z = f.blocks{3};
    
    % Make sure all components are real-valued.
    if ( ~isreal(x) || ~isreal(y) || ~isreal(z) )
        error('CHEBFUN:CHEBMATRIX:compose:complex3', ...
            'Can compose only with a CHEBMATRIX with real-valued CHEBFUN entries.')
    end
    
    if ( isa(op, 'chebfun3') )
        f = chebfun(@(t) op(feval(x, t), feval(y, t), feval(z, t)), x.domain);
        
    elseif ( isa(op, 'chebfun3v') )
        % Call compose for each component of OP.
        F = compose(f, op.components{1});
        for jj = 2:op.nComponents
            F = [ F; compose(f, op(jj)) ]; % op.components{jj})
        end
        f = F;
        
    else
        error('CHEBFUN:CHEBMATRIX:compose:entries3', ...
            'Can compose a 3 by 1 CHEBMATRIX only with a CHEBFUN3 or CHEBFUN3V.')
    end
    
else
    error('CHEBFUN:CHEBMATRX:compose:rows', ...
        'Can compose only with an 2 by 1 or 3 by 1 chebmatrix.')
    
end

end
