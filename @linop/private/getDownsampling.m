function [d,dRow,dVar] = getDownsampling(L)

[m,n] = size(L);
dRow = getRowDiffOrders(L);
dVar = getColDiffOrders(L);
% Each variable requires a downsampling contribution equal to its
% differential order. The total diff. orders of the variables might
% not equal the total of the rows.
totalRow = sum( dRow(~isnan(dRow)) );
totalVar = sum( dVar(~isnan(dVar)) );
if (totalRow == totalVar)
    d = dRow;
else
    % The only reasonable thing is to spread out the D.O. reductions
    % as evenly as possible among the rows.
    d = NaN(1,m);
    isOp = ~isnan(dRow);
    numOp = sum(isOp);
    k = floor( totalVar / numOp );
    d(isOp) = k;  % even distribution of order
    i = rem(totalVar,numOp);  % leftover to be made up
    if i > 0
        % Increment the first i rows
        incr = find(isOp,i);
        d(incr) = d(incr) + 1;
    end
end

end
