function f = compose(f, op, varargin)
% COMPOSE   Compose command for CHEBFUN3 objects. 
%   F = COMPOSE(F, OP) returns a CHEBFUN3 that approximates OP(F).
% 
%   F = COMPOSE(F, OP, G) returns a CHEBFUN3 that approximates OP(F, G).
%   This command is a wrapper for the CHEBFUN3 constructor. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ( nargin == 2 ) && ( nargin(op) == 1 ) )
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
