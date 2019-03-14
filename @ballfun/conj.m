function g = conj(f)
%CONJ  Complex conjugate of a BALLFUN.
%   CONJ(F) is the complex conjugate of F.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

g = compose( f, @conj ); 

end
