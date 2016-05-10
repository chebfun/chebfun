function f = mrdivide( f, g )
% /   Right scalar divide for SEPARABLEAPPROX objects.
%
%   F/C divides the SEPARABLEAPPROX F by a scalar C.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(g, 'double') )
    f.pivotValues = f.pivotValues * g;
else
    error('CHEBFUN:SEPARABLEAPPROX:mrdivide:mrdivide', ...
        'Not supported. Did you mean ./ ?');
end

end
