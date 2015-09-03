function h = times(f, g)
% .*    Pointwise multiplication for SEPARABLEAPPROX objects.
%
%   F.*G multiplies SEPARABLEAPPROX objects F and G. Alternatively F or G could be a
%   double.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(f, 'separableApprox') )    % SEPARABLEAPPROX .* ???
    
    if ( isa(g, 'double') )  % SEPARABLEAPPROX .* DOUBLE
        h = mtimes(f, g);
    elseif ( isa( g, 'separableApprox') )
        bol = domainCheck(f, g);
        if ( bol )
            h = compose( f, @times, g); 
        else
            error('CHEBFUN:SEPARABLEAPPROX:times:domain', 'Inconsistent domains');
        end
    else
        error('CHEBFUN:SEPARABLEAPPROX:times:unknown', ...
            ['Undefined function ''mtimes'' for input arguments of type %s ' ...
            'and %s.'], class(f), class(g));
    end
    
else
    h = times(g, f);
end

end
