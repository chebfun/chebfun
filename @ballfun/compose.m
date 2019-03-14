function f = compose(f, op, g )
%COMPOSE   Composition for BALLFUN.
%   COMPOSE(F, OP) returns a BALLFUN representing OP(F) where F is also a
%   BALLFUN object, and OP is a function handle.
%
%   COMPOSE(F, OP, G) returns a BALLFUN representing OP(F, G) where F and G
%   are BALLFUN objects, and OP is a function handle.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(op) )
    return
elseif ( isempty(f) )
    return
else
    % A very simple compose function for now: 
    if ( nargin == 2 )
        op = @(r, l, th) feval(op, feval(f, r, l, th, 'spherical'));
    else
        if ( isempty(g) )
            op = f;
        else
            op = @(r, l, th) feval(op, feval(f, r, l, th, 'spherical'), feval(g, r, l, th, 'spherical'));
        end
    end    

    % Make BALLFUN object:
    f = ballfun( op, 'spherical' );
end
end
    