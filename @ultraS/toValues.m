function fx = toValues(disc, f)
%TOVALUES  Convert a CHEBFUN to its ULTRAS discretization.

% Input should be a single (perhaps piecewise smooth) chebfun. 

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

dom = disc.domain;
numInts = numel(dom) - 1;
dim = disc.dimension;

f = restrict(f, dom);

c = cell(numInts, 1);
for k = 1:numInts
    c{k} = flipud(chebpoly(f.funs{k}, dim(k)));
end
c = cell2mat(c);
S = convert(disc, 0, disc.outputSpace);
fx = S*c;

end
