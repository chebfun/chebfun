function f = compose(f, op, varargin)
%COMPOSE   Compose command for CHEBFUN3 objects.
%   F = COMPOSE(F, OP) returns a CHEBFUN3 that approximates OP(F).
%
%   F = COMPOSE(F, OP, G) returns a CHEBFUN3 that approximates OP(F, G).
%   This command is a wrapper for the CHEBFUN3 constructor.
%
%   F = COMPOSE(F, G) with a CHEBFUN G with one column returns a CHEBFUN3 that
%   approximates G(F).  If G has 2 or 3 columns, the result is a CHEBFUN3V.
%
%   F = COMPOSE(F, G) for a CHEBFUN2 or CHEBFUN2V G returns G(F) interpreted as
%   G(real(F), imag(F)), regardless whether F is real or complex valued.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(op, 'chebfun') )
    % Composition OP(f) of CHEBFUN3 object f and CHEBFUN OP
    
    if ( ~isreal(f) )
        error('CHEBFUN:CHEBFUN3:COMPOSE:Complex', ...
            'Composition of a CHEBFUN and a complex CHEBFUN3 is not defined.')
    end
    
    if ( length(op.domain) > 2 )
        % If OP has several pieces, OP(CHEBFUN3) might be inaccurate.
        warning('CHEBFUN:CHEBFUN3:compose:pieces', ...
            ['The composition of a CHEBFUN with several pieces and a CHEBFUN3\n', ...
            'might be inaccurate.']);
    end
    
    % Check that image(f) is contained in domain(OP).
    vals = minandmax3est(f);        % Estimate of image(f).
    tol = 100 * chebfun3eps * max(vscale(f), vscale(op)) * ...
        norm(f.domain, inf);        % Tolerance.
    if ( ~isSubset(vals, op.domain, tol) )
        error('CHEBFUN:CHEBFUN3:COMPOSE:DomainMismatch', ...
            'OP(F) is not defined, since image(F) is not contained in domain(OP).')
    end
    
    nColumns = size(op, 2);
    if ( nColumns <= 3 )
        % If f is periodic, then OP(f) should be periodic:
        pref = chebfunpref;
        if ( isPeriodicTech(f) )
            pref.tech = get(f.rows(:,1).funs{1}, 'tech');
        end
        
        % Compute first entry:
        opcolumn = op(:,1);
        F = chebfun3(@(x,y,z) opcolumn(feval(f, x, y, z)), f.domain, pref);
        % Add additional components:
        for jj = 2:nColumns
            opcolumn = op(:,jj);
            F = [ F; chebfun3(@(x,y,z) opcolumn(feval(f, x, y, z)), ...
                f.domain, pref) ];
        end
        f = F;
        
    else
        % The CHEBFUN object OP has a wrong number of columns.
        error('CHEBFUN:CHEBFUN3:COMPOSE:Columns', ...
            'The CHEBFUN object must have 1, 2, or 3 columns.')
        
    end
    
elseif ( isa(op, 'chebfun2') || isa(op, 'chebfun2v') )
    % Interpret OP(f) as OP(real(f), imag(f))
    F = [real(f); imag(f)];
    f = compose(F, op);
    
elseif ( ( nargin == 2 ) && ( nargin(op) == 1 ) )
    % OP has one input variable.
    
    % If f is periodic, then OP(f) should be periodic:
    pref = chebfunpref;
    if ( isPeriodicTech(f) )
        pref.tech = get(f.rows(:,1).funs{1}, 'tech');
    end
    
    % Call constructor:
    f = chebfun3(@(x,y,z) op(feval(f, x, y, z)), f.domain, pref);
    
elseif ( ( nargin == 3 ) && ( nargin(op) == 2 ) )
    % OP has two input variables.
    
    g = varargin{1};
    if ( isa(g, 'double') )     % promote
        g = chebfun3(g, f.domain);
    end
    
    if ( isa(f, 'double') )     % promote
        f = chebfun3(f, g.domain);
    end
    
    % Call constructor:
    f = chebfun3(@(x,y,z) op(feval(f, x, y, z), feval(g, x, y, z)), ...
        f.domain);
    
else
    % Not sure what to do, error:
    error('CHEBFUN:CHEBFUN3:COMPOSE:OP', 'NARGIN(OP) not correct.')
    
end

end
