function d = getEqnDiffOrders(L)
[m,n] = size(L);
d = zeros(1,m);
for i = 1:m
    for j = 1:n
        block = L.blocks{i,j};
        if isa(block,'linopOperator')
            d(i) = max(d(i),block.diffOrder);
        elseif isa(block,'linopFunctional') || isnumeric(block)
            d(i) = NaN;
        end
    end
end

end
