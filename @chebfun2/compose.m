function f = compose(f, op, varargin)
%COMPOSE   Compose command for CHEBFUN2 objects.
%   F = COMPOSE(F, OP) returns the CHEBFUN2 that approximates OP(F).
%
%   F = COMPOSE(F, OP, G) returns the CHEBFUN2 that approximates OP(F,G).
%   This command is a wrapper for the CHEBFUN2 constructor.
%
%   F = COMPOSE(F, G) with a CHEBFUN G with one column returns a CHEBFUN2 that
%   approximates G(F).  If G has 2 or 3 columns, the result is a CHEBFUN2V.
%
%   F = COMPOSE(F, G) for a CHEBFUN2 or CHEBFUN2V G returns G(F) interpreted as
%   G(real(F), imag(F)), regardless whether F is real or complex valued.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(op, 'chebfun') )
    % Composition OP(f) of CHEBFUN2 object f and CHEBFUN OP
    
    if ( ~isreal(f) )
        error('CHEBFUN:CHEBFUN2:COMPOSE:Complex', ...
            'Composition of a CHEBFUN and a complex CHEBFUN2 is not defined.')
    end
    
    if ( length(op.domain) > 2 )
        % If OP has several pices, OP(CHEBFUN2) might be inaccurate.
        warning('CHEBFUN:CHEBFUN2:compose:pieces', ...
            ['The composition of a CHEBFUN with several pieces and a CHEBFUN2\n', ...
            'might be inaccurate.']);
    end
    
    % TO DO: Should we add the following domain check?
    %     % Check that image(f) is contained in domain(OP).
    %     vals = minandmax2(f);
    %     if ( ( vals(1) < op.domain(1) ) || ( vals(2) > op.domain(2) ) )
    %         error('CHEBFUN:CHEBFUN2:COMPOSE:DomainMismatch', ...
    %             'Composition op(f) not defined, because the image of f is not contained in the domain of op.')
    %     end
    
    nColumns = size(op, 2);
    if ( nColumns <= 3 )
        % Compute first entry:
        opcolumn = op(:,1);
        
        if ( isPeriodicTech(f) )
            % OP(f) should be periodic if f is:
            F = chebfun2(@(x,y) opcolumn(feval(f, x, y)), f.domain, 'trig');
            % Add additional components:
            for jj = 2:nColumns
                opcolumn = op(:,jj);
                F = [ F; chebfun2(@(x,y) opcolumn(feval(f, x, y)), ...
                    f.domain, 'trig') ];
            end
            
        else
            F = chebfun2(@(x,y) opcolumn(feval(f, x, y)), f.domain);
            % Add additional components:
            for jj = 2:nColumns
                opcolumn = op(:,jj);
                F = [ F; chebfun2(@(x,y) opcolumn(feval(f, x, y)), f.domain) ];
            end
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
    
    % Call constructor:
    if ( isPeriodicTech(f) )
        f = chebfun2(@(x,y) op(feval(f, x, y)), f.domain, 'trig');
    else
        f = chebfun2(@(x,y) op(feval(f, x, y)), f.domain);
    end
    
elseif ( ( nargin == 3 ) && ( nargin(op) == 2 ) )
    % OP has two input variables.
    
    g = varargin{1};
    if ( isa(g, 'double') )     % promote
        g = chebfun2(g, f.domain);
    end
    
    if ( isa(f, 'double') )     % promote
        f = chebfun2(f, g.domain);
    end
    
    % Call constructor:
    f = chebfun2(@(x,y) op(feval(f, x, y), feval(g, x, y)), f.domain);
    
else
    % Not sure what to do, error:
    error('CHEBFUN:CHEBFUN2:COMPOSE:OP', 'NARGIN(OP) not correct.')
    
end

end