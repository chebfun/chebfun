function F = real(F)
%REAL   Real part of a CHEBFUN3V object.
%   REAL(F) returns the real part of a CHEBFUN3V object.
%
% See also CHEBFUN3V/IMAG and CHEBFUN3V/CONJ.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(F) )
    return
end

% Take real part of each component:
F.components = cellfun(@real, F.components, 'UniformOutput', false);

end