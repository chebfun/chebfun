function [B,P] = resize(disc,A,m,n)

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
if ( m == n )
    B = A;
    P = eye(size(A,1));
    return
end
        

domain = disc.domain;
numint = disc.numIntervals;
P = cell(1,numint);

for k = 1:numint
    disc.domain = domain(k:k+1);
    disc.dimension = m(k);
    xOut = points(disc,1);
    disc.dimension = n(k);
    xIn = points(disc,2);
    P{k} = barymat(xOut,xIn);
end

P = blkdiag(P{:});
B = P*A;

end
