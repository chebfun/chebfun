function f = compose(f, op, varargin)
%COMPOSE   Compose command for SPHEREFUN objects.
%   F = COMPOSE(F, OP)  returns the SPHEREFUN that approximates OP(F).
% 
%   F = COMPOSE(F, OP, G)  returns the SPHEREFUN that approximates OP(F).
%
%   F = COMPOSE(F, G) with a CHEBFUN G with one column returns a SPHEREFUN that
%   approximates G(F).  If G has 3 columns, the result is a SPHEREFUNV. If G is
%   a CHEBFUN2 or CHEBFUN2V, the composition is interpreted as G(real(F),
%   imag(F)).
%
%   This command is a wrapper for the SPHEREFUN constructor.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(op) )
    return
elseif ( isempty(f) )
    f = op;
    
elseif ( isa(op, 'chebfun') )
    % Composition OP(f) of SPHEREFUN object f and CHEBFUN OP
    
    if ( length(op.domain) > 2 )
        % If OP has several pieces, OP(SPHEREFUN) might be inaccurate.
        warning('CHEBFUN:SPHEREFUN:compose:pieces', ...
            ['The composition of a CHEBFUN with several pieces and a SPHEREFUN\n', ...
            'might be inaccurate.']);
    end
    
    % Check that image(f) is contained in domain(OP).
    vals = minandmax2est(f);        % Estimate of image(f).
    tol = 100 * chebfun2eps * max(vscale(f), vscale(op)) * ...
            norm(f.domain, inf);    % Tolerance.
    if ( ~isSubset(vals, op.domain, tol) )
        error('CHEBFUN:SPHEREFUN:COMPOSE:DomainMismatch', ...
            'OP(F) is not defined, since image(F) is not contained in domain(OP).')
    end
    
    nColumns = size(op, 2);
    if ( nColumns == 1 )
        % Call constructor:
        f = spherefun(@(x,y) op(feval(f, x, y)), f.domain);
        
    elseif ( nColumns == 3 )
        % Extract columns of the CHEBFUN OP:
        op1 = op(:,1);
        op2 = op(:,2);
        op3 = op(:,3);
        
        % Call constructor:
        f = spherefunv(@(x,y) op1(feval(f, x, y)), ...
            @(x,y) op2(feval(f, x, y)), @(x,y) op3(feval(f, x, y)));
        
    else
        % The CHEBFUN object OP has a wrong number of columns.
        error('CHEBFUN:SPHEREFUN:COMPOSE:Columns', ...
            'The CHEBFUN object must have 1 or 3 columns.')
        
    end
    
elseif ( isa(op, 'chebfun2') )
    % Composition OP(f) of SPHEREFUN object f and CHEBFUN2 OP, interpreted as
    % OP(real(f), imag(f)).  For now SPHEREFUNS are real, but this might change
    % in the future.
        
    % Check that image(f) is contained in domain(OP).
    vals = minandmax2est(f);        % Estimate of image(f).
    tol = 100 * chebfun2eps * max(vscale(f), vscale(op)) * ...
            norm(f.domain, inf);    % Tolerance.
    if ( ~isSubset(vals, op.domain(1:2), tol) )
        error('CHEBFUN:SPHEREFUN:COMPOSE:DomainMismatch2', ...
            'OP(F) is not defined, since image(F) is not contained in domain(OP).')
    end
    
    % Call constructor:
    f = spherefun(@(x,y) op(feval(real(f), x, y), feval(imag(f), x, y)), ...
        f.domain);
    
elseif ( isa(op, 'chebfun2v') )
    % Composition OP(f) of SPHEREFUN object f and CHEBFUN2V OP with three
    % components, interpreted as OP(real(f), imag(f)).
    % For now SPHEREFUNS are real, but this might change in the future.
        
    % Check that OP has three components:
    if ( op.nComponents ~= 3 )
        error('CHEBFUN:SPHEREFUN:COMPOSE:C2VnComponents', ...
            'The Chebfun2v objects must have three components.')
    end
    
    % Get the components:
    op1 = op(1);
    op2 = op(2);
    op3 = op(3);
    
    % Check that image(f) is contained in domain(OP).
    vals = minandmax2est(f);        % Estimate of image(f).
    tol = 100 * chebfun2eps * max(vscale(f), vscale(op)) * ...
            norm(f.domain, inf);    % Tolerance.
    if ( ~isSubset(vals, op1.domain(1:2), tol) )
        error('CHEBFUN:SPHEREFUN:COMPOSE:DomainMismatch2v', ...
            'OP(F) is not defined, since image(F) is not contained in domain(OP).')
    end
    
    % Get real and imaginary parts of f:
    realf = real(f);
    imagf = imag(f);
    
    % Call constructor:
    f = spherefunv(@(x,y) op1(feval(realf, x, y), feval(imagf, x, y)), ...
        @(x,y) op2(feval(realf, x, y), feval(imagf, x, y)), ...
        @(x,y) op3(feval(realf, x, y), feval(imagf, x, y)));
    
elseif ( nargin == 2 && nargin(op) == 1 )
    % OP has one input variable.
    
    % Call constructor: 
    f = spherefun(@(x,y) op( feval(f, x, y) ), f.domain);
    
elseif ( nargin == 3 && nargin(op) == 2 )
    % OP has two input variables. 
    
    g = varargin{1}; 
    if ( isa(g, 'double') )     % promote
        g = spherefun(@(x,y,z) g + 0*x, f.domain);
    end
    
    if ( isa(f, 'double') )     % promote
        f = spherefun(@(x,y,z) f + 0*x, g.domain); 
    end
    
    % Call constructor: 
    f = spherefun(@(x,y) op( feval(f, x, y), feval(g, x, y) ), f.domain);
    
else
    % Not sure what to do, error: 
    error('CHEBFUN:SPHEREFUN:COMPOSE:OP', 'NARGIN(OP) not correct.') 
end

end 