function [L, S] = quasi2diffmat(disc)
%QUASI2USDIFFMAT(DISC)   Convert DISC.coeffs to a differential operator.
%   L = QUASI2USDIFFMAT(DISC) returns a matrix L corresponding to the
%   differential operator C{1}*D^[m] + ... C{m+1}*I.
%
% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get info from DISC:
c = fliplr(disc.coeffs);
if ( isempty(disc.outputSpace) )
    disc.outputSpace = size(c, 2) - 1;
end
dim = disc.dimension;

L = sparse(sum(dim), sum(dim));
for j = 1:size(c, 2)
    L = L + mult(disc, c{j}) * diff(disc, j - 1);
end

if ( nargout > 1 )
    S = [];
end

end
