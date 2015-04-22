function h = times(f, g)
% .*    Pointwise multiplication for LOWRANKAPPROX objects.
%
%   F.*G multiplies LOWRANKAPPROX objects F and G. Alternatively F or G could be a
%   double.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(f, 'lowrankapprox') )    % LOWRANKAPPROX .* ???
    
    if ( isa(g, 'double') )  % LOWRANKAPPROX .* DOUBLE
        h = mtimes(f, g);
    elseif ( isa( g, 'lowrankapprox') )
        bol = domainCheck(f, g);
        if ( bol )
            h = compose( f, @times, g); 
        else
            error('CHEBFUN:LOWRANKAPPROX:times:domain', 'Inconsistent domains');
        end
    else
        error('CHEBFUN:LOWRANKAPPROX:times:unknown', ...
            ['Undefined function ''mtimes'' for input arguments of type %s ' ...
            'and %s.'], class(f), class(g));
    end
    
else
    h = times(g, f);
end

end
