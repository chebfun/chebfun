function h = times(f, g)
% .*   CHEBFUN3 multiplication.
%
%   F.*G multiplies CHEBFUN3 objects F and G. Alternatively F or G could 
%   be a double.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(f, 'chebfun3') )     % CHEBFUN3 .* ???
    
    if ( isa(g, 'double') )  % CHEBFUN3 .* DOUBLE
        h = mtimes(f, g);
    elseif ( isa( g, 'chebfun3') )
        bool = domainCheck(f, g);
        if ( bool )
            h = chebfun3(@(x, y, z) feval(f, x, y, z) .* feval(g, x, y, z), ...
                f.domain, 'fiberDim', 3);
        else
           error('CHEBFUN:CHEBFUN3:times:domain', 'Inconsistent domains');
       end
    else
        error('CHEBFUN:CHEBFUN3:times:unknown', ...
            ['Undefined function ''mtimes'' for input arguments of type %s ' ...
            'and %s.'], class(f), class(g));
    end

else 
   h = times(g, f); 
end