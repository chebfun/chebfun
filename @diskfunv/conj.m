function F = conj( F )
%CONJ Complex conjugate of a DISKFUNV.
%   CONJ(F) returns the complex conjugate of F. For a complex F, CONJ(F) =
%   REAL(F) - i*IMAG(F).
%
%   Since only real DISKFUNV objects are supported, the result of this
%   function is always the original function F.
%
%   See also DISKFUNV/IMAG, DISKFUNV/REAL.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty( F ) )
    return
end

% Take conj part of each component:
F.components = cellfun( @conj, F.components, 'UniformOutput', false );

end
