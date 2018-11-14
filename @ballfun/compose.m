function f = compose(f, op, g )
%COMPOSE   Composition for BALLFUN.
%   COMPOSE(F, OP) returns a BALLFUN representing OP(F) where F is also a
%   BALLFUN object, and OP is a function handle.
%
%   COMPOSE(F, OP, G) returns a BALLFUN representing OP(F, G) where F and G
%   are BALLFUN objects, and OP is a function handle.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% A very simple compose function for now: 
if ( nargin == 2 )
    op = @(r, l, th) feval(op, feval(f, r, l, th));
else
    op = @(r, l, th) feval(op, feval(f, r, l, th), feval(g, r, l, th));
end    

% Make BALLFUN object:
f = ballfun( op, 'polar' );
end
    