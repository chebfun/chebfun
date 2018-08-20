function g = conj(f)
% CONJ Conjugate of a BALLFUN function
%   CONJ(f) is the conjugate of the BALLFUN function f

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

g = real(f)-1i*imag(f);
end
