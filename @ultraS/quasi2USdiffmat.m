function L = quasi2USdiffmat(disc)

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.


c = disc.coeffs;
dom = disc.domain;
dim = disc.dimension;
outputSpace = disc.outputSpace;

if isempty(outputSpace)
    diffOrder = size(c, 2) - 1;
    outputSpace = diffOrder;
end

dummy = ultraS([]);
dummy.domain = dom;
dummy.dimension = dim;
c = fliplr(c);

L = 0*speye(sum(dim));
for j = 1:size(c, 2)
    %form D^(j-1) term.
    L = L + convert(dummy, j-1, outputSpace)*mult(dummy, c{j}, j-1)*diff(dummy, j - 1);
end


end