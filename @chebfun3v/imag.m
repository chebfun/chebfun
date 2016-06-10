function F = imag( F )
%IMAG   Imaginary part of a CHEBFUN3V.
%
%   IMAG(F) returns the imaginary part of a CHEBFUN3V.
%
%   See also CONJ, REAL.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty(F) )
    return
end

% Take imag part of each component:
F.components = cellfun(@imag, F.components, 'UniformOutput', false);

end