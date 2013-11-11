clc
clear classes

%% Building blocks
dom = [-2 2];
I = linop.eye(dom);
D = linop.diff(dom);
Z = linop.zeros(dom);
x = chebfun('x', dom);
u = chebfun('x.^2', dom);
U = linop.mult(u);   

%%
Dexact = [ 
      -2.750000000000000   3.414213562373095  -1.000000000000000   0.585786437626905  -0.250000000000000
  -0.853553390593274   0.353553390593274   0.707106781186548  -0.353553390593274   0.146446609406726
   0.250000000000000  -0.707106781186548                   0   0.707106781186548  -0.250000000000000
  -0.146446609406726   0.353553390593274  -0.707106781186548  -0.353553390593274   0.853553390593274
   0.250000000000000  -0.585786437626905   1.000000000000000  -3.414213562373095   2.750000000000000
];

%% Collocation discretizations
zero1 = norm( discretize(I,5) - eye(5) )
zero2 = norm( discretize(D,5) - Dexact );
xx = colloc2.points(5,dom);
zero3 = norm( discretize(U,5) - diag(u(xx)) )

%% Building blocks
dom = [-2 1 1.5 2];
I = linop.eye(dom);
D = linop.diff(dom);
Z = linop.zeros(dom);
x = chebfun('x', dom);
u = chebfun('x.^2', dom);
U = linop.mult(u);   

%% Collocation discretizations
zero4 = norm( discretize(I,5) - eye(15) )
xx = blockColloc2.points([5 5 5],dom);
zero5 = norm( discretize(U,5) - diag(u(xx)) )
