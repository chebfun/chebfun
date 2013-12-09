function G = cheb2quasi(F)

if ( numel(F) > 1 )
    G = F;
    return
end

F = num2cell(F);

numCols = numel(F);
G(numCols) = chebfun();
for k = 1:numCols
    G(k) = F{k};
end

end