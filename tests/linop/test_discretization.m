function pass = test_discretization

%% Building blocks
dom = [-2 2];
I = chebmatrix( operatorBlock.eye(dom) );
D = chebmatrix( operatorBlock.diff(dom) );
u = chebfun('x.^2', dom);
U = chebmatrix( operatorBlock.mult(u) );   

%%
Dexact = [ 
      -2.750000000000000   3.414213562373095  -1.000000000000000   0.585786437626905  -0.250000000000000
  -0.853553390593274   0.353553390593274   0.707106781186548  -0.353553390593274   0.146446609406726
   0.250000000000000  -0.707106781186548                   0   0.707106781186548  -0.250000000000000
  -0.146446609406726   0.353553390593274  -0.707106781186548  -0.353553390593274   0.853553390593274
   0.250000000000000  -0.585786437626905   1.000000000000000  -3.414213562373095   2.750000000000000
];

%% Collocation discretizations
err(1) = norm( matrix(I,5) - eye(5) );
err(2) = norm( matrix(D,5) - Dexact );
xx = chebpts(5, dom);
err(3) = norm( matrix(U,5) - diag(u(xx)) );

%% Building blocks
dom = [-2 1 1.5 2];
I = chebmatrix( operatorBlock.eye(dom) );
u = chebfun('x.^2', dom);
U = chebmatrix( operatorBlock.mult(u) );  
n = [5 5 5];

%% Collocation discretizations
err(4) = norm( matrix(I,n) - eye(sum(n)) );
xx = chebpts(n,dom);
err(5) = norm( matrix(U,n) - diag(u(xx)) );

pass = err < 1e-9;

end
