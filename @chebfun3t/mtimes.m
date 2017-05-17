function h = mtimes(f, g)
%*  Scalar multiplication for CHEBFUN3T objects.
%   c*F or F*c multiplies a CHEBFUN3T F by a scalar c.
%
% See also CHEBFUN3T/TIMES.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(f, 'chebfun3t') )           % CHEBFUN3T * ???
    
    if ( isa(g, 'double') )         % CHEBFUN3T * DOUBLE
        if ( numel(g) == 1 )
            if ( g == 0 ) 
                h = chebfun3t(@(x,y,z) 0, f.domain);
            else
                h = f;
                h.coeffs = h.coeffs .* g;
                h.vscale = h.vscale .* abs(g);
            end
        else
            error('CHEBFUN:CHEBFUN3T:mtimes:size', 'Sizes are inconsistent.');
        end
        
    else
        error('CHEBFUN:CHEBFUN3T:mtimes:unknown', ...
            ['Undefined function ''mtimes'' for input arguments of type %s ' ...
            'and %s.'], class(f), class(g));
        
    end
    
elseif ( isa(f, 'double') && isa(g, 'chebfun3t') )  % DOUBLE * CHEBFUN3T
    h = mtimes(g.', f.').';
else
    error('CHEBFUN:CHEBFUN3T:mtimes:size', 'Sizes are inconsistent.');
end

end