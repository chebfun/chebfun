function discjk = extractBlock(disc, j, k)
%EXTRACTBLOCK   Extract the j-k block from a discretization.
%   DISCJK = EXTRACTBLOCK(DISC, J, K) extracts information relating to the
%   J-K block of the dscretization DISC and returns a new discretization
%   DISCJK. DISCJK.dimension will be modified by the amount specified in
%   DISC.dimAdjust(j,k).

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

discjk = extractBlock@opDiscretization(disc, j, k);

% Extract the coefficients and output space:
discjk.coeffs = disc.coeffs{j, k};
discjk.outputSpace = disc.outputSpace(j);

end
