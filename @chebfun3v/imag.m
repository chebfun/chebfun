function F = imag(F)
%IMAG   Imaginary part of a CHEBFUN3V object.
%   IMAG(F) returns the imaginary part of a CHEBFUN3V.
%
% See also CHEBFUN3V/CONJ and CHEBFUN3V/REAL.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty(F) )
    return
end

% Take imaginary part of each component:
F.components = cellfun(@imag, F.components, 'UniformOutput', false);

end