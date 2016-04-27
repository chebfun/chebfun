function F = real( F )
%REAL  Real part of a SPHEREFUNV.
%   REAL(F) returns the SPHEREFUNV representing the real part.
%
%   Since only real SPHEREFUNV objects are supported, the result of this
%   function is always just the original function F.
%
%   See also SPHEREFUNV/CONJ, SPHEREFUNV/IMAG.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( F ) )
    return
end

% Take real part of each component:
F.components = cellfun( @real, F.components, 'UniformOutput', false );

end
