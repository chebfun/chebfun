function h = times(f, g)
%.*   CHEBFUN3T pointwise multiplication.
%   F.*G multiplies CHEBFUN3T objects F and G. Alternatively F or G could be
%   a double.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(f, 'chebfun3t') )     % CHEBFUN3T .* ???
    
    if ( isa(g, 'double') )    % CHEBFUN3T .* DOUBLE
        h = mtimes(f, g);
    elseif ( isa(g, 'chebfun3t') )
        sameDom = domainCheck(f, g);
        if ( sameDom )
            h = chebfun3t(@(x, y, z) feval(f, x, y, z).*feval(g, x, y, z), ...
                f.domain);
        else
           error('CHEBFUN:CHEBFUN3T:times:domain', 'Inconsistent domains');
       end
    else
        error('CHEBFUN:CHEBFUN3T:times:unknown', ...
            ['Undefined function ''times'' for input arguments of type %s ' ...
            'and %s.'], class(f), class(g));
    end

else 
   h = times(g, f); 
end

end