function S = summat(N)
%SUMMAT   Definite integral matrix for TRIGSPEC.
%   S = SUMMAT(N) returns the definite integral NxN matrix for the TRIGSPEC
%   class.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

Z = zeros(N, 1);
Z(floor(N/2) + 1) = 1;
S = spdiags(Z, 0, N, N);

end