function h = mldivide(f, g)
%\   Left divide for CHEBFUN3T objects.
%   MLDIVIDE(F, G) computes left dividison for CHEBFUN3T objects F and G. 
%   Only allowed to divide by scalars.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% If either f or g are empty then return an empty CHEBFUN3 object.
if ( isempty(f) )
    h = f;
    return
    
elseif ( isempty(g) )
    h = g;
    return 
end

if ( ~isa(f, 'double') )
    error('CHEBFUN:CHEBFUN3T:mldivide:nonScalar', ...
        'Division must be by a scalar.');
end

% Left divide:
h = ldivide(f, g);

end