function D = diff(A,m)

d = A.domain;
n = dim(A);

if m == 0
    D = eye(sum(n));
else
    numIntervals = length(d)-1;
    
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
