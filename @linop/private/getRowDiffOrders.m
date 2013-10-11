function d = getRowDiffOrders(L)

[m,n] = size(L);
d = zeros(1,m);

for i = 1:m
    for j = 1:n
        block = L.operator.blocks{i,j};
        if isa(block,'operatorBlock')
            d(i) = max(d(i),block.diffOrder);
        elseif isa(block,'functionalBlock') || isnumeric(block)
            d(i) = NaN;
        end
    end
end

end
