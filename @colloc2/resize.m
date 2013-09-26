function B = resize(A,m,n,domain)
numint = length(m);
P = cell(1,numint);
for k = 1:numint
    xOut = colloc2.points(m(k),domain(k:k+1));
    xIn = colloc2.points(n(k),domain(k:k+1));
    P{k} = barymat(xOut,xIn);
end
B = blkdiag(P{:})*A;
end
