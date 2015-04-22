function f = mrdivide( f, g )
% /   Right scalar divide for LOWRANKAPPROX objects.
%
%   F/C divides the LOWRANKAPPROX F by a scalar C.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(g, 'double') )
    f.pivotValues = f.pivotValues * g;
else
    error('CHEBFUN:LOWRANKAPPROX:mrdivide:mrdivide', ...
        'Not supported. Did you mean ./ ?');
end

end
