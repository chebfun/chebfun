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
        g = diskfun(g, f.domain); 
    end
    
    if ( isa( f, 'double' ) )     % promote
        f = diskfun(f, g.domain); 
    end
    
    % Call constructor: 
    f = diskfun(@(x,y) op( feval(f, x, y,1), feval(g, x, y,1) ), f.domain); %polar
else
    % Not sure what to do, error: 
    error('CHEBFUN:DISKFUN:COMPOSE:OP', 'NARGIN(OP) not correct.')
    
end

end 