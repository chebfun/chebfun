function F = ctranspose(F)
%'   Complex conjugate transpose of a CHEBFUN.
%   F' is the complex conjugate transpose of the CHEBFUN F.
%
% See also TRANSPOSE.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

for k = 1:numel(F)
    F(k) = transpose(conj(F(k)));
end

end
