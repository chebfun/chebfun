function f = ctranspose( f )
%'	 Complex conjugate transpose of a CHEBFUN2.
%   F' is the complex conjugate transpose of F.
%   G = CTRANSPOSE(F) is called for the syntax F'.  
%
% See also CONJ, TRANSPOSE.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Take conj part:
f = conj( f ); 

% Call transpose:
f = transpose( f ); 

end
