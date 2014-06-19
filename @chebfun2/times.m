function h = times(f, g)
% .*   CHEBFUN2 multiplication.
%   F.*G multiplies CHEBFUN2 objects F and G. Alternatively F or G could be a
%   double.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(f, 'chebfun2') )    % CHEBFUN2 .* ???
    
    if ( isa(g, 'double') )  % CHEBFUN2 .* DOUBLE
        h = mtimes(f, g);
    elseif ( isa( g, 'chebfun2') )
        bol = domainCheck(f, g);
        if ( bol )
            h = chebfun2(@(x, y) feval(f, x, y).*feval(g, x, y), f.domain);
        else
            error('CHEBFUN:CHEBFUN2:times:domain', 'Inconsistent domains');
        end
    else
        error('CHEBFUN:CHEBFUN2:times:unknown', ...
            ['Undefined function ''mtimes'' for input arguments of type %s ' ...
            'and %s.'], class(f), class(g));
    end
    
else
    h = times(g, f);
end

end
