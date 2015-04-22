function h = mldivide(f, g)
%\      Left divide for LOWRANKAPPROX objects.
%
% Left divide for a LOWRANKAPPROX. Only allowed to divide by scalars.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty( f ) || isempty( g ) )
    h = lowrankapprox();
    return 
end

if ( ~isa(f, 'double') )
    error('CHEBFUN:LOWRANKAPPROX:mldivide:nonScalar', ...
        'Division must be by a scalar.');
end

% Left divide:
h = ldivide( f, g );

end
