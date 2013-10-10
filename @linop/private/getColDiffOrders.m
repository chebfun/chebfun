function d = getColDiffOrders(L)

[m,n] = size(L);
d = zeros(1,n);

for j = 1:n
    for i = 1:m
        block = L.operator.blocks{i,j};
        if isa(block,'operatorBlock')
            d(j) = max(d(j),block.diffOrder);
        elseif isa(block,'functionalBlock') || isnumeric(block)
            d(i) = NaN;
        end
    end
end

end
