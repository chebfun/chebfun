function varargout = eig(N,k)
%EIG

if nargin < 2 
    k = 6; 
end


rect = N.domain;
m = 25; n = 25;
f = chebfun2(0, rect);

CC = constructDiscretisation(N,f,m,n,1); % Construct discretisation.

% make massive n^2 by n^2 matrix.
sz = (size(CC{1,1},1))*(size(CC{2,1},2));

A = spalloc(sz, sz, m*sz + sz);  P = diag((1:n).^(-2));
for jj = 1:size(CC,1)
    A1 = CC{jj,2}; A1([end-1,end],:) = [ones(1,n);(-1).^(0:n-1)];
    B1 = CC{jj,1}; B1([end-1,end],:) = [ones(1,n);(-1).^(0:n-1)];
    A = A + kron(A1,B1);
end

S1 = spconvermat(m,0,N.xorder); S1([end-1,end],:) = zeros(2,n);
S2 = spconvermat(m,0,N.yorder); S2([end-1,end],:) = zeros(2,n);
S = kron(S1,S2);

% P = diag((1:m*n).^(-2));
if nargout > 1 
    [VV,X] = eig(full(A), full(S));
    X = diag(X); 
else
    X = sort(eig(full(A), full(S)));
end

idx = find(~isinf(X) & abs(X)<100 & imag(X)<sqrt(eps));
X = real(X(idx));

[ignored,I] = sort(abs(X));
X = X(I); 
X = X(1:min(length(X),k));

V = cell(1,length(X));
if nargout > 1
    VV = VV(:,idx); 
    VV = VV(:,I);
    VV  = VV(:,1:length(X));
    for jj = 1:length(X)
        mat = reshape(VV(:,jj),size(CC{1,1},1),size(CC{2,1},2));
        mat = chebfun2(polyval(rot90(mat,2)), rect);
        V{1,jj} = mat; 
    end
    varargout = {V,X};  
else
   varargout = {X};  
end
        
end

    
function c = polyval(c)
% Convert bivariate chebyshev coeffs to values on a Chebyshev grid.
[sy, sx] = size(c);
for j = 1:sx
    c(:,j) = chebpolyval(c(:,j));
end
for k = 1:sy
    c(k,:) = chebpolyval(c(k,:).').';
end
end