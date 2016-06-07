function f = compose(f, op, varargin)
%COMPOSE    compose command for SPHEREFUN objects.
% 
%   F = COMPOSE(F, OP)  returns the SPHEREFUN that approximates OP(F).
% 
%   F = COMPOSE(F, OP, G)  returns the SPHEREFUN that approximates OP(F).
%
%   This command is a wrapper for the SPHEREFUN constructor.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


if ( isempty(op) )
    return
elseif ( isempty(f) )
    f = op;
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