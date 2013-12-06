function S = sum(disc)

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
n = disc.dimension;
d = disc.domain;

numint = disc.numIntervals;

S = cell(1,numint);
for k = 1:numint
    S{k} = chebtech2.quadwts(n(k)) * (d(k+1)-d(k))/2;
end
S = cat(2,S{:});

end
