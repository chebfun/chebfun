function f = real(f)
%REAL   Real part of a SINGFUN.
%   REAL(F) is the real part of F.
%
% See also IMAG, ISREAL, CONJ.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

% Compute the real part of the smooth part of F.
f.smoothPart = real(f.smoothPart);

% If F is imaginary, then the smooth part will be zero due to the line above.
% In this case remove redundant singularities from what is being returned.
if ( iszero(f.smoothPart) )
    f.exponents = [0, 0];    
    % NH: f = f.smoothPart;
end

end