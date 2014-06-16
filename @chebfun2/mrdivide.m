function f = mrdivide( f, g )
% /   Right scalar divide
%   F/C divides the CHEBFUN2 F by a scalar C.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(g, 'double') )
    f.pivotValues = f.pivotValues * g;
else
    error('CHEBFUN:CHEBFUN2:mrdivide:mrdivide', ...
        'Not supported. Did you mean ./ ?');
end

end
