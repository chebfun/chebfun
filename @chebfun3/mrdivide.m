function f = mrdivide(f, g)
%/   Right scalar divide for CHEBFUN3 objects.
%
%    F/C divides the CHEBFUN3 object F by a scalar C.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(g, 'double') )
    f.core = f.core / g;
else
    error('CHEBFUN:CHEBFUN3:mrdivide:mrdivide', ...
        'Not supported. Did you mean ./ ?');
end

end