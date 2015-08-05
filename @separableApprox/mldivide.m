function h = mldivide(f, g)
%\      Left divide for SEPARABLEAPPROX objects.
%
% Left divide for a SEPARABLEAPPROX. Only allowed to divide by scalars.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% If either f or g are empty then return an empty SEPARABLEAPPROX object.
if ( isempty( f ) )
    h = f;
    return;
elseif ( isempty( g ) )
    h = g;
    return 
end

if ( ~isa(f, 'double') )
    error('CHEBFUN:SEPARABLEAPPROX:mldivide:nonScalar', ...
        'Division must be by a scalar.');
end

% Left divide:
h = ldivide( f, g );

end
