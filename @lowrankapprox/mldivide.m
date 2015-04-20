function h = mldivide(f, g)
%\   CHEBFUN2 left divide.
%
% Left divide for a CHEBFUN2. Only allowed to divide by scalars.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) || isempty(g) )
    h = chebfun2;
    return
end

if ( ~isa(f, 'double') )
    error('CHEBFUN:CHEBFUN2:mldivide:nonScalar', ...
        'Division must be by a scalar.');
end

% Left divide:
h = ldivide( f, g );

end
