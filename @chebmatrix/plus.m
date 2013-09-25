function C = plus(A,B)

[m,n] = size(A);
C = cell(m,n);

for i = 1:m
    for j = 1:n
        C{i,j} = A.blocks{i,j} + B.blocks{i,j};
    end
end

C = chebmatrix(C);

end
