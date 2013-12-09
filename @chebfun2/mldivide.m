function h = mldivide(f,g)
%/   Chebfun2 left divide.
%
% Left divide for a chebfun2. Only allowed to divide by scalars.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if isempty(f) || isempty(g)
    h = chebfun2;
    return;
end

if ~isa(f,'double')
    error('CHEBFUN2:MRDIVIDE:NONSCALAR','Division must be by a scalar.');
end


% Left divide.
h = ldivide(f,g);

end