function B = resize(disc,A,m,n)

domain = disc.domain;
numint = disc.numIntervals;
P = cell(1,numint);

for k = 1:numint
    disc.domain = domain(k:k+1);
    disc.dimension = m(k);
    xOut = points(disc);
    disc.dimension = n(k);
    xIn = points(disc);
    P{k} = barymat(xOut,xIn);
end

B = blkdiag(P{:})*A;

end
