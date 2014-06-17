function out = isreal(f)
%ISREAL   True for real DELTAFUN.
%   ISREAL(F) returns TRUE if the FUNPART and all the delta functions of F 
%   are real and FALSE otherwise.
%
% See also REAL, IMAG.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Check if the funPart is real:
if ( ~isreal(f.funPart) )
    out = 0;
    return
end

% Check the delta functions:
out = isreal(f.deltaMag);

end
