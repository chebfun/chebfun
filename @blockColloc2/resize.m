function B = resize(A,m,n,domain, difforder)

numint = length(m);
P = cell(1,numint);

for k = 1:numint
    xOut = blockColloc2.points(m(k),domain(k:k+1));
    xIn = blockColloc2.points(n(k),domain(k:k+1));
    P{k} = barymat(xOut,xIn);
end

B = blkdiag(P{:})*A;

end
