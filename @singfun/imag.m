function f = imag(f)
%IMAG   Imaginary part of a SINGFUN.
%   IMAG(F) is the imaginary part of F.
%
% See also REAL, ISREAL, CONJ.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check for empty arguments:
if ( isempty(f) )
    return;
end

% Compute the imaginary part of the smooth part:
f.smoothPart = imag(f.smoothPart);

% Return a SMOOTHFUN object if F is smooth:
if ( issmooth(f) )
    f = f.smoothPart;
end
    
end
