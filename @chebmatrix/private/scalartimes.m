function C = scalartimes(A,z)
[m,n] = size(A);
C = cell(m,n);
for i = 1:m
    for j = 1:n
        C{i,j} = z*A.blocks{i,j};
    end
end
C = chebmatrix(C);
end
