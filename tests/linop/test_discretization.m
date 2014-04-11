function pass = test_discretization

%% Building blocks
dom = [-2 2];
I = chebmatrix( operatorBlock.eye(dom) );
D = chebmatrix( operatorBlock.diff(dom) );
u = chebfun('x.^2', dom);
U = chebmatrix( operatorBlock.mult(u) );   

x5 = chebpts(5, dom, 2);
x6 = chebpts(6, dom, 2);
y5 = chebpts(5, dom, 1);
P55 = barymat(y5, x5);
P56 = barymat(y5, x6);
D56 = .5*P56*diffmat(6);

%% Collocation discretizations
err(1) = norm( matrix(I,5) - P55 );
err(2) = norm( matrix(D,5) - D56 );
xx = chebpts(5, dom);
err(3) = norm( matrix(U,5) - P55*diag(u(xx)) );

%% Building blocks
dom = [-2 1 1.5 2];
I = chebmatrix( operatorBlock.eye(dom) );
u = chebfun('x.^2', dom);
U = chebmatrix( operatorBlock.mult(u) );  
n = [5 5 5];

x = chebpts(n, dom, 2);
y = chebpts(n, dom, 1);
P55 = cell(3,1);
for k = 1:3
    idx = (k-1)*5+(1:5);
    P55{k} = barymat(y(idx), x(idx));
end
P55 = blkdiag(P55{:});

%% Collocation discretizations
err(4) = norm( matrix(I,n) - P55 );
xx = chebpts(n, dom);
err(5) = norm( matrix(U,n) - P55*diag(u(xx)) );

pass = err < 1e-9;

end
