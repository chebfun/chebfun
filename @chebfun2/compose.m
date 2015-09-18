function f = compose( f, op, varargin )
% COMPOSE     compose command for CHEBFUN2 objects. 
% 
%  F = COMPOSE(F, OP )  returns the CHEBFUN2 that approximates OP(F).
% 
%  F = COMPOSE(F, OP, G )  returns the CHEBFUN2 that approximates OP(F,G).
%
% This command is a wrapper for the CHEBFUN2 constructor. 

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 2 && nargin(op) == 1)
    % OP has one input variable.
    
    % Call constructor: 
    f = chebfun2(@(x,y) op( feval(f, x, y) ), f.domain);
    
elseif ( nargin == 3 && nargin(op) == 2 )
    % OP has two input variables. 
    
    g = varargin{1}; 
    if ( isa( g, 'double' ) )     % promote
        g = chebfun2(g, f.domain); 
    end
    
    if ( isa( f, 'double' ) )     % promote
        f = chebfun2(f, g.domain); 
    end
    
    % Call constructor: 
    f = chebfun2(@(x,y) op( feval(f, x, y), feval(g, x, y) ), f.domain);
    
else
    % Not sure what to do, error: 
    error('CHEBFUN:CHEBFUN2:COMPOSE:OP', 'NARGIN(OP) not correct.')
    
end

end 