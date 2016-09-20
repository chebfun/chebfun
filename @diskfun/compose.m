function f = compose( f, op, varargin )
% COMPOSE   compose command for DISKFUN objects. 
% 
%   F = COMPOSE(F, OP ) returns the DISKFUN that approximates OP(F).
% 
%   F = COMPOSE(F, OP, G ) returns the DISKFUN that approximates OP(F, G).
%
% This command is a wrapper for the DISKFUN constructor. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 2 && nargin(op) == 1)
    % OP has one input variable.
    f = diskfun(@(x,y) op( feval(f, x, y, 'polar') ), 'polar');
elseif ( nargin == 3 && nargin(op) == 2 )
    % OP has two input variables. 
    g = varargin{1}; 
    if ( isa( g, 'double' ) )     % promote
        g = diskfun(@(x,y) g + 0*x); 
    end
    
    if ( isa( f, 'double' ) )     % promote
        f = diskfun(@(x,y) f + 0*x); 
    end
    
    f = diskfun(@(x,y) op( feval(f, x, y, 'polar'), feval(g, x, y, 'polar') ), 'polar'); 
else
    % Not sure what to do, error: 
    error('CHEBFUN:DISKFUN:COMPOSE:OP', 'NARGIN(OP) not correct.')
end

end 

