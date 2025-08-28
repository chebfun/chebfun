function f = compose(f, op, varargin)
%COMPOSE   Compose command for CHEBFUN2 objects.
%   F = COMPOSE(F, OP) returns the CHEBFUN2 that approximates OP(F).
%
%   F = COMPOSE(F, OP, G) returns the CHEBFUN2 that approximates OP(F,G).
%   This command is a wrapper for the CHEBFUN2 constructor.
%
%   F = COMPOSE(F, G) with a CHEBFUN G with one column returns a CHEBFUN2 
%   that approximates G(F). If G has 2 or 3 columns, the result is a CHEBFUN2V.
%
%   F = COMPOSE(F, G) for a CHEBFUN2 or CHEBFUN2V G returns G(F) 
%   interpreted as G(real(F), imag(F)), regardless of whether F is real or 
%   complex valued.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(op, 'chebfun') )
    % Composition OP(f) of CHEBFUN2 object f and CHEBFUN OP
    
    if ( ~isreal(f) )
        error('CHEBFUN:CHEBFUN2:COMPOSE:Complex', ...
            'Composition of a CHEBFUN and a complex CHEBFUN2 is not defined.')
    end
    
    if ( length(op.domain) > 2 )
        % If OP has several pieces, OP(CHEBFUN2) might be inaccurate.
        warning('CHEBFUN:CHEBFUN2:compose:pieces', ...
            ['The composition of a CHEBFUN with several pieces and a CHEBFUN2\n', ...
            'might be inaccurate.']);
    end
    
    % Check that image(f) is contained in domain(OP):
    vals = minandmax2est(f);        % Estimate of image(f).
    tol = 100 * chebfun2eps * max(vscale(f), vscale(op)) * ...
        norm(f.domain, inf);        % Tolerance.
    if ( ~isSubset(vals, op.domain, tol) )
        error('CHEBFUN:CHEBFUN2:COMPOSE:DomainMismatch', ...
            'OP(F) is not defined, since image(F) is not contained in domain(OP).')
    end
    
    nColumns = size(op, 2);
    if ( nColumns <= 3 )
        % If f is periodic in x or y, then OP(f) should inherit that:
        funtype = 'cheb';
        isTrig = [isPeriodicTech(f.rows), isPeriodicTech(f.cols)];
        if ( isTrig(1) && isTrig(2) )
            funtype = 'trig';
        elseif ( isTrig(1) && ~isTrig(2) )
            funtype = 'trigx';
        elseif ( ~isTrig(1) && isTrig(2) )
            funtype = 'trigy';
        end
        
        % Compute first entry:
        opcolumn = op(:,1);
        F = chebfun2(@(x,y) opcolumn(feval(f, x, y)), f.domain, funtype);
        % Add additional components:
        for jj = 2:nColumns
            opcolumn = op(:,jj);
            F = [F; chebfun2(@(x,y) opcolumn(feval(f, x, y)), f.domain, funtype)];
        end
        f = F;
        
    else
        % The CHEBFUN object OP has a wrong number of columns.
        error('CHEBFUN:CHEBFUN2:COMPOSE:Columns', ...
            'The CHEBFUN object must have 1, 2, or 3 columns.')
    end
    
elseif ( isa(op, 'chebfun2') || isa(op, 'chebfun2v') )
    % Interpret OP(f) as OP(real(f), imag(f))
    F = [real(f); imag(f)];
    f = compose(F, op);
    
elseif ( ( nargin == 2 ) && ( nargin(op) == 1 ) )
    % OP has one input variable.
    
    % If f is periodic in x or y, then OP(f) should inherit that:
    funtype = 'cheb';
    isTrig = [isPeriodicTech(f.rows), isPeriodicTech(f.cols)];
    if ( isTrig(1) && isTrig(2) )
        funtype = 'trig';
    elseif ( isTrig(1) && ~isTrig(2) )
        funtype = 'trigx';
    elseif ( ~isTrig(1) && isTrig(2) )
        funtype = 'trigy';
    end
    
    % Call constructor:
    f = chebfun2(@(x,y) op(feval(f, x, y)), f.domain, funtype);
    
elseif ( ( nargin == 3 ) && ( nargin(op) == 2 ) )
    % OP has two input variables.
    
    g = varargin{1};

    % OP(f,g) should be periodic in x or y if both f and g are periodic in
    % that direction. Note that a scalar is periodic in both directions. 
    funtype = 'cheb';
    isTrigf = [false, false]; %whether f is trig in x and y
    isTrigg = [false, false]; %whether g is trig in x and y
    if ( isa(f, 'chebfun2') )
        isTrigf = [isPeriodicTech(f.rows), isPeriodicTech(f.cols)];
    elseif ( isscalar(f) )
        isTrigf = [true, true];
    end
    if ( isa(g, 'chebfun2') )
        isTrigg = [isPeriodicTech(g.rows), isPeriodicTech(g.cols)];
    elseif ( isscalar(g) )
        isTrigg = [true, true];
    end
    if ( isTrigf(1) && isTrigg(1) && isTrigf(2) && isTrigg(2))
        funtype = 'trig';
    elseif ( isTrigf(1) && isTrigg(1) )
        funtype = 'trigx';
    elseif ( isTrigf(2) && isTrigg(2) )
        funtype = 'trigy';
    end

    if ( isa(g, 'double') )     % promote
        g = chebfun2(g, f.domain);
    end
    
    if ( isa(f, 'double') )     % promote
        f = chebfun2(f, g.domain);
    end
    
    %Domain check
    if ( ~domainCheck(f, g) )
        error('CHEBFUN:CHEBFUN2:COMPOSE:DomainIncompatible', ...
            'The CHEBFUN2 objects do not share the same domain')
    end
    
    % Call constructor:
    f = chebfun2(@(x,y) op(feval(f, x, y), feval(g, x, y)), f.domain, funtype);
    
else
    % Not sure what to do, error:
    error('CHEBFUN:CHEBFUN2:COMPOSE:OP', 'NARGIN(OP) not correct.')   
end

end
