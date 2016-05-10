function F = real( F )
%REAL  Real part of a CHEBFUN2V.
%   REAL(F) returns the CHEBFUN2V representing the real part.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( F ) )
    return
end

% Take real part of each component:
F.components = cellfun( @real, F.components, 'UniformOutput', false );

end
