function f = compose( f, op, varargin )
% COMPOSE     compose command for DISKFUN objects. 
% 
%  F = COMPOSE(F, OP )  returns the DISKFUN that approximates OP(F).
% 
%  F = COMPOSE(F, OP, G )  returns the DISKFUN that approximates OP(F).
%
% This command is a wrapper for the DISKFUN constructor. 

if ( nargin == 2 && nargin(op) == 1)
    % OP has one input variable.
    
    % Call constructor: 
    f = diskfun(@(x,y) op( feval(f, x, y) ), f.domain);
    
elseif ( nargin == 3 && nargin(op) == 2 )
    % OP has two input variables. 
    
    g = varargin{1}; 
    if ( isa( g, 'double' ) )     % promote
        g = diskfun(@(x,y) g + 0*x, f.domain); 
    end
    
    if ( isa( f, 'double' ) )     % promote
        f = diskfun(@(x,y) f + 0*x, g.domain); 
    end
    
    % Call constructor: 
    f = diskfun(@(x,y) op( feval(f, x, y), feval(g, x, y) ), f.domain); %polar
else
    % Not sure what to do, error: 
    error('CHEBFUN:DISKFUN:COMPOSE:OP', 'NARGIN(OP) not correct.')
    
end

end 

