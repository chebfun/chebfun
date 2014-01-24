function F = real( F )
%REAL  real part of a chebfun2v.
%
% REAL(F) returns the chebfun2v representing the real part.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Empty check: 
if ( isempty( F ) )
    return
end

% Take real part of each component:
F.components = cellfun( @real, F.components, 'UniformOutput', false );

end