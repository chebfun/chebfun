function fx = toValues(disc, f)
%TOVALUES   Convert a CHEBFUN to its ULTRAS discretization.
%   C = TOVALUES(DISC, F) converts the (perhaps piecewise smooth) chebfun F
%   to coefficients C for use by an ULTRAS discretization DISC.

% Copyright 20134 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

dom = disc.domain;
numInts = disc.numIntervals();
dim = disc.dimension;

% Restrict F to DOM (e.g. if we need to introduce more breakpoints in F)
f = restrict(f, dom);

% Loop through the pieces of F and obtain their Chebyshev coefficients, stored
% in the cell-array C.
c = cell(numInts, 1);
for k = 1:numInts
    c{k} = flipud(chebpoly(f.funs{k}, dim(k)));
end
c = cell2mat(c); 

% Create a conversion operator for the ultrasphericial method, using the
% appropriate discretisation and output space.
S = convert(disc, 0, disc.outputSpace);

% Apply the conversion operator to the vector C, so that we get coefficients in
% the correct order Chebyshev basis.
fx = S*c;

end
