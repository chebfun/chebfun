function h = mldivide(f, g)
%\      Left divide for CHEBFUN3 objects.
%
% Left divide for a CHEBFUN3. Only allowed to divide by scalars.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% If either f or g are empty then return an empty CHEBFUN3 object.
if ( isempty(f) )
    h = f;
    return;
    
elseif ( isempty(g) )
    h = g;
    return 
end

if ( ~isa(f, 'double') )
    error('CHEBFUN:CHEBFUN3:mldivide:nonScalar', ...
        'Division must be by a scalar.');
end

% Left divide:
h = ldivide(f, g);

end