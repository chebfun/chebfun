function D = diff(A, m)
%DIFF    Differentiation operator for ultraspherical spectral method
%   D = DIFF(A, m) returns a differentiation operator for the ultrapsherical
%   spectral method that represents the mth derivative.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Useful information:
d = A.domain;
n = A.dimension;

if ( m == 0 )          
    % 0th order derivative is easy!
    D = speye(sum(n));
    
else
    numIntervals = length(d) - 1;
    
    % Find the diagonal blocks:
    blocks = cell(numIntervals);
    
    % Construct the blocks:
    for k = 1:numIntervals
        len = d(k+1) - d(k);
        blocks{k} = ultraS.diffmat(n(k), m) * (2/len)^m;
    end
    
    % Assemble:
    D = blkdiag(blocks{:});
end

end
