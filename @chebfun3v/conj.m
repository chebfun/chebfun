function F = conj(F)
%CONJ   Complex conjugate of a CHEBFUN3V.
%   CONJ(F) returns the complex conjugate of F. For a complex F, we have
%   CONJ(F) = REAL(F) - i*IMAG(F).
%
% See also CHEBFUN3V/REAL and CHEBFUN3V/IMAG.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty(F) )
    return
end

% Take conj part of each component:
F.components = cellfun(@conj, F.components, 'UniformOutput', false);

end