function h = times(f, g)
% .*   CHEBFUN2 multiplication.
%   F.*G multiplies CHEBFUN2 objects F and G. Alternatively F or G could be a
%   double.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( isa(f, 'chebfun2') )    % CHEBFUN2 .* ???
    
    if ( isa(g, 'double') )  % CHEBFUN2 .* DOUBLE
        h = mtimes(f, g);
    elseif ( isa( g, 'chebfun2') )
        bol = domainCheck(f, g);
        if ( bol )
            h = chebfun2(@(x, y) feval(f, x, y).*feval(g, x, y), f.domain);
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