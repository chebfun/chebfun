function f = real(f)
%REAL   Real part of a SINGFUN.
%   REAL(F) is the real part of F.
%
% See also IMAG, ISREAL, CONJ.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check for empty arguments:
if ( isempty(f) )
    return;
end

% Compute the real part of the smooth part of F.
f.smoothPart = real(f.smoothPart);

% Return a SMOOTHFUN object if F is smooth:
if ( issmooth(f) )
    f = f.smoothPart;
end

end
