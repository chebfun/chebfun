function F = ctranspose(F)
%'    Complex conjugate transpose.
%   F' is the complex conjugate transpose of the CHEBFUN F.
%
% See also TRANSPOSE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

F = transpose(conj(F));

end
