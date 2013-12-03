function D = diff(A, m)
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
d = A.domain;
n = A.dimension;
if ( m == 0 )
    D = speye(sum(n));
else
    numIntervals = length(d) - 1;
    % Find the diagonal blocks.
    blocks = cell(numIntervals);
    for k = 1:numIntervals
        len = d(k+1) - d(k);
        blocks{k} = diffmat(n(k), m) * (2/len)^m;
    end
    % Assemble.
    D = blkdiag(blocks{:});
end
end
