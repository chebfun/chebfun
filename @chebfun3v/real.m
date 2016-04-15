function F = real( F )
%REAL  Real part of a CHEBFUN3V.
%   REAL(F) returns the CHEBFUN3V representing the real part.
%   See also CONJ, IMAG.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( F ) )
    return
end

% Take real part of each component:
F.components = cellfun( @real, F.components, 'UniformOutput', false );

end