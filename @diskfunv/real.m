function F = real( F )
%REAL  Real part of a DISKFUNV.
%   REAL(F) returns the DISKFUNV representing the real part of F.
%
%   Since all DISKFUNV objects are real-valued, the result of this
%   function is always the original function F.
%
% See also DISKFUNV/CONJ, DISKFUNV/IMAG.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( F ) )
    return
end

% Take real part of each component:
F.components = cellfun( @real, F.components, 'UniformOutput', false );

end