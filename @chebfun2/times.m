function h = times(f, g)
% .*   CHEBFUN2 times.
%

if ( isa(f, 'chebfun2') )    % CHEBFUN2 .* ???
    if ( isa(g, 'double') )  % CHEBFUN2 .* DOUBLE
        h = mtimes(f, g);
    elseif ( isa( g, 'chebfun2') )
        [bol, dom] = domain_check(f, g);
        if ( bol )
            h = chebfun2(@(x, y) feval(f, x, y).*feval(g, x, y), dom);
        else
            error('CHEBFUN2:TIMES:DOMAIN', 'Inconsistent domains');
        end
    else
        error('CHEBFUN:MTIMES:unknown', ...
            ['Undefined function ''mtimes'' for input arguments of type %s ' ...
            'and %s.'], class(f), class(g));
    end
else
    h = times(g, f);
end

end

function [bol, domain] = domain_check(f, g)
% Check that the domains of f and g are the same.

Fdom = f.domain;
Gdom = g.domain;
Fscl = max(abs(Fdom));
Gscl = max(abs(Gdom));

bol = ( norm(Fdom - Gdom) < eps * max(Fscl, Gscl) );
domain = (Fdom + Gdom)/2;

end