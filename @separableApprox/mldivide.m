function h = mldivide(f, g)
%\      Left divide for SEPARABLEAPPROX objects.
%
% Left divide for a SEPARABLEAPPROX. Only allowed to divide by scalars.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty( f ) || isempty( g ) )
    h = separableApprox();
    return 
end

if ( ~isa(f, 'double') )
    error('CHEBFUN:SEPARABLEAPPROX:mldivide:nonScalar', ...
        'Division must be by a scalar.');
end

% Left divide:
h = ldivide( f, g );

end
