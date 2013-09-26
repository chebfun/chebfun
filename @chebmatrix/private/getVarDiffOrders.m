function d = getVarDiffOrders(L)
[m,n] = size(L);
d = zeros(1,n);
for j = 1:n
    for i = 1:m
        block = L.blocks{i,j};
        if isa(block,'linopOperator')
            d(j) = max(d(j),block.diffOrder);
        elseif isa(block,'linopFunctional') || isnumeric(block)
            d(i) = NaN;
        end
    end
end

end
