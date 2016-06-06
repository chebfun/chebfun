function f = mrdivide(f, g)
%/   Right scalar divide for CHEBFUN3 objects.
%   MRDIVIDE(F, C) or F/C divides the CHEBFUN3 object F by a scalar C.
%
% See also CHEBFUN3/MLDIVIDE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% If either f or g are empty then return an empty CHEBFUN3 object.
if ( ( isempty(f) ) ||  ( isempty(g) ) )
    return
end

if ( isa(g, 'double') )
    if ( isscalar(g) )
        if ( g == 0 )
            error('CHEBFUN:CHEBFUN3:mrdivide:divisionByZero', ...
                'Division by zero.')
        end
    end
    
    f.core = f.core / g;
else
    error('CHEBFUN:CHEBFUN3:mrdivide:mrdivide', ...
        'Not supported. Did you mean ./ ?');
end

end