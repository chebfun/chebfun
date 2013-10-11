function C = cumsum(A,m)

d = A.domain;
n = dim(A);

if m == 0
    C = eye(sum(n));
else
    numIntervals = length(d)-1;
    
    % Find the diagonal blocks.
    blocks = cell(numIntervals);
    for k = 1:numIntervals
        len = d(k+1) - d(k);
        blocks{k} = cumsummat(n(k)) * (len/2);
    end
    
    % Assemble.
    C = blkdiag(blocks{:});
    
    % Each subinterval also contributes to the integrals in all the
    % subintervals to its right, creating a triangular structure.
    offset = 0;
    for k = 1:numIntervals
        % Grab the weights for the integral using all of this
        % subinterval.
        row = offset + n(k);
        cols = offset + (1:n(k));
        last = C(row,cols);
        % Copy it to add to the ones that follow.
        offset = row;
        C(offset+1:end,cols) = repmat(last,[sum(n)-offset,1]);
    end
    C = C^m;
end
end
