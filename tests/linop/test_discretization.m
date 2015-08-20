function pass = test_discretization

% TODO: Tests 2,3,5 assume a chebcolloc2 discretization.

pref = cheboppref();
pref.discretization = @chebcolloc2;

%% Building blocks
dom = [-2 2];
I = chebmatrix( operatorBlock.eye(dom) );
D = chebmatrix( operatorBlock.diff(dom) );
u = chebfun('x.^2', dom);
U = chebmatrix( operatorBlock.mult(u) );   

D55 = [  
  -2.750000000000000   3.414213562373094  -1.000000000000000   0.585786437626905  -0.250000000000000
  -0.853553390593274   0.353553390593274   0.707106781186548  -0.353553390593274   0.146446609406726
   0.250000000000000  -0.707106781186548                   0   0.707106781186548  -0.250000000000000
  -0.146446609406726   0.353553390593274  -0.707106781186548  -0.353553390593274   0.853553390593274
   0.250000000000000  -0.585786437626905   1.000000000000000  -3.414213562373094   2.750000000000000];

%% Collocation discretizations
err(1) = norm( matrix(I,5,pref) - eye(5) );
err(2) = norm( matrix(D,5,pref) - D55 );
xx = chebpts(5, dom);
err(3) = norm( matrix(U,5,pref) - diag(u(xx)) );

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

%% Collocation discretizations
err(4) = norm( matrix(I,n,pref) - eye(sum(n)) );
xx = chebpts(n, dom);
err(5) = norm( matrix(U,n,pref) - diag(u(xx)) );
pass = err < 1e-9;

end
