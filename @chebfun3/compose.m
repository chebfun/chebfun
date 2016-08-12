function f = compose(f, op, varargin)
%COMPOSE   Compose command for CHEBFUN3 objects.
%   F = COMPOSE(F, OP) returns a CHEBFUN3 that approximates OP(F).
%
%   F = COMPOSE(F, OP, G) returns a CHEBFUN3 that approximates OP(F, G).
%   This command is a wrapper for the CHEBFUN3 constructor.
%
%   F = COMPOSE(F, G) with a CHEBFUN G returns a CHEBFUN3 object that
%   approximates G(F).  If G has 2 or 3 columns, the result is a CHEBFUN3V.
%
%   F = COMPOSE(F, G) with a CHEBMATRIX G of size n by 1 (n = 1,2,3) returns a
%   CHEBFUN3 or CHEBFUN3V object that approximates G(F).

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(op, 'chebfun') )
    % Composition OP(f) of CHEBFUN3 object f and CHEBFUN OP
    
    if ( ~isreal(f) )
        error('CHEBFUN:CHEBFUN3:COMPOSE:Complex', ...
            'Composition of a CHEBFUN and a complex CHEBFUN3 is not defined.')
    end
    
    % TO DO: Should we add the following domain check?
    %     % Check that image(f) is contained in domain(OP).
    %     vals = minandmax3(f);
    %     if ( ( vals(1) < op.domain(1) ) || ( vals(2) > op.domain(2) ) )
    %         error('CHEBFUN:CHEBFUN3:COMPOSE:DomainMismatch', ...
    %             'Composition op(f) not defined, because the image of f is not contained in the domain of op.')
    %     end
    
    nColumns = size(op, 2);
    if ( nColumns == 1 )
        % Call constructor:
        f = chebfun3(@(x,y,z) op(feval(f, x, y, z)), f.domain);
        
    elseif ( nColumns == 2 );
        % Extract columns of OP:
        op1 = op(:,1);
        op2 = op(:,2);
        
        % Call constructor:
        f = chebfun3v(@(x,y,z) op1(feval(f, x, y, z)), ...
            @(x,y,z) op2(feval(f, x, y, z)), f.domain);
        
    elseif ( nColumns == 3 )
        % Extract columns of OP:
        op1 = op(:,1);
        op2 = op(:,2);
        op3 = op(:,3);
        
        % Call constructor:
        f = chebfun3v(@(x,y,z) op1(feval(f, x, y, z)), ...
            @(x,y,z) op2(feval(f, x, y, z)), ...
            @(x,y,z) op3(feval(f, x, y, z)), f.domain);
        
    else
        % The CHEBFUN object OP has a wrong number of columns.
        error('CHEBFUN:CHEBFUN3:COMPOSE:Columns', ...
            'The CHEBFUN object must have 1, 2, or 3 columns.')
        
    end
    
elseif ( isa(op, 'chebmatrix') )
    % Composition OP(f) of the CHEBFUN3 object f and the CHEBMATRIX OP.
    
    if ( ~isreal(f) )
        error('CHEBFUN:CHEBFUN3:COMPOSE:Complex', ...
            'Composition of a CHEBMATRIX and a complex CHEBFUN3 is not defined.')
    end
    
    [m, n] = size(op);
    if ( ( m > 3 ) || ( n > 1 ) )
        error('CHEBFUN:CHEBFUN3:COMPOSE:ChebmatrixFormat', ...
            'The CHEBMATRIX must by n by 1 with n=1,2,3.')
        
    else
        F = f;
        % Compose each entry of OP with F.
        f = compose(F, op{1});
        for iRow = 2:m
            f = [f; compose(F, op{iRow})];
        end
        
    end
    
elseif ( isa(op, 'chebfun2') || isa(op, 'chebfun2v') )
    % Interpret OP(f) as OP(real(f), imag(f))
    F = [real(f); imag(f)];
    f = compose(F, op);
    
elseif ( ( nargin == 2 ) && ( nargin(op) == 1 ) )
    % OP has one input variable.
    
    % Call constructor:
    f = chebfun3(@(x,y,z) op(feval(f, x, y, z)), f.domain);
    
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
