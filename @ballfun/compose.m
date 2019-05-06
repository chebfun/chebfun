function f = compose(f, op, g )
%COMPOSE   Compose command for BALLFUN objects.
%   COMPOSE(F, OP) returns a BALLFUN that approximates OP(F).
%
%   COMPOSE(F, OP, G) returns a BALLFUN approximates OP(F, G).

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
        % Make BALLFUN object:
        f = ballfun( op, 'spherical' );
    else
        if isempty(g)
            f = g;
        else
            op = @(r, l, th) feval(op, feval(f, r, l, th, 'spherical'), feval(g, r, l, th, 'spherical'));
            % Make BALLFUN object:
            f = ballfun( op, 'spherical' );
        end
    end    
end
end
    