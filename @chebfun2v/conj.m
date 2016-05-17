function F = conj( F )
%CONJ Complex conjugate of a CHEBFUN2V.
%   CONJ(F) returns the complex conjugate of F. For a complex F, CONJ(F) =
%   REAL(F) - i*IMAG(F).
%
% See also REAL, IMAG. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty( F ) )
    return
end

% Take conj part of each component:
F.components = cellfun( @conj, F.components, 'UniformOutput', false );

end
