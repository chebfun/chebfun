function f = real(f)
%REAL	Real part of a CLASSICFUN.
%   REAL(F) is the real part of F.
%
% See also IMAG, ISREAL, CONJ.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Assign the ONEFUN to be the real part of the ONEFUN of the input:
f.onefun = real(f.onefun);

end
