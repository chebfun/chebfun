function D = diff(disc,m)

d = disc.domain;
n = disc.dimension;

if m == 0
    D = eye(sum(n));
else
    numIntervals = disc.numIntervals;
    
    % Find the diagonal blocks.
    blocks = cell(numIntervals);
    for k = 1:numIntervals
        len = d(k+1) - d(k);
        blocks{k} = diffmat(n(k),m) * (2/len)^m;
    end
    
    % Assemble.
    D = blkdiag(blocks{:});
end
end
