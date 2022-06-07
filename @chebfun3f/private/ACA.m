%% Adaptive Cross Approximation with full pivoting
function [Ac, At, Ar, rowInd, colInd] = ACA(A, tol, maxIter)
%A \approx Ac inv(At) Ar'
Ac = [];
Ar = [];
At = [];
rowInd = [];
colInd = [];
Aoriginal = A;

for iter = 1:maxIter
    
    [error,I2] = max(abs(A(:)));
    if isempty(error) || error < tol
        Ac = Aoriginal(:,colInd);
        Ar = Aoriginal(rowInd,:)';
        At = Aoriginal(rowInd,colInd);
        return
    end
    
    [I,J] = ind2sub(size(A), I2);
    rowInd = [rowInd, I];
    colInd = [colInd, J];
    
    A = A-A(:,J)*A(I,:)./A(I,J);
end
Ac = Aoriginal(:,colInd);
Ar = Aoriginal(rowInd,:)';
At = Aoriginal(rowInd,colInd);
end