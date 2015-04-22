function f = compose( f, op, varargin )
% COMPOSE     compose command for CHEBFUN2 objects. 
% 
%  F = COMPOSE(F, OP )  returns the CHEBFUN2 that approximates OP(F).
% 

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 2 )
    
    f = chebfun2(@(x,y) op( feval(f, x, y) ), f.domain);
    
elseif ( nargin == 3 && nargin(op) == 2 )
    
    g = varargin{1}; 
    if ( isa( g, 'double' ) )     %promote
        g = chebfun2(g, f.domain); 
    end
    
    if ( isa( f, 'double' ) )    %promote
        f = chebfun2(f, g.domain); 
    end
    
    f = chebfun2(@(x,y) op( feval(f, x, y), feval(g, x, y) ), f.domain);
    
end

end 