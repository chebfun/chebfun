function f = mrdivide( f, g )
% /   Right scalar divide
%   F/C divides the CHEBFUN2 F by a scalar C.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( isa(g, 'double') )
    f.pivotValues = f.pivotValues * g;
else
    error('CHEBFUN2:mrdivide', 'Not supported. Did you mean ./ ?');
end

end
