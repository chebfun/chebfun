function idx = getDealiasingIndexes(~, N, nVars)
%GETGDEALIASINGINDEXES  % Returns indexes for dealiasing procedure (2/3-rule).
%
% See also SPINOPSPHERE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

toOne = floor(N/2) + 1 - ceil(N/6):floor(N/2) + ceil(N/6);
idx = false(N, N);
idx(toOne, toOne) = 1;
idx = repmat(idx, nVars, 1);

end