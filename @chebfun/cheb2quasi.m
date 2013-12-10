function G = cheb2quasi(F)

if ( numel(F) > 1 )
    G = F;
    return
end

F = num2cell(F);

numCols = numel(F);
for k = numCols:-1:1
    G(k) = F{k};
end

end