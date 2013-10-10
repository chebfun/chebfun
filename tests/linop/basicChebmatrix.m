clc
clear classes

%% Building blocks
dom = [-2 -0.5 1 2];
I = linop.eye(dom);
D = linop.diff(dom);
Z = linop.zeros(dom);
x = chebfun('x', dom);
u = sin(x.^2);
U = linop.diag(u);   

D5 = [ 
  -5.499999999999999   6.828427124746189  -2.000000000000000   1.171572875253810  -0.500000000000000
  -1.707106781186547   0.707106781186547   1.414213562373095  -0.707106781186548   0.292893218813452
   0.500000000000000  -1.414213562373095                   0   1.414213562373095  -0.500000000000000
  -0.292893218813452   0.707106781186548  -1.414213562373095  -0.707106781186547   1.707106781186547
   0.500000000000000  -1.171572875253810   2.000000000000000  -6.828427124746189   5.499999999999999
];

%%
A = chebmatrix( {I,Z;D,U} );
M = discretize(A,5);
DD = blkdiag(2/1.5*D5,2/1.5*D5,2/1*D5);
[xx,ww] = blockColloc2.points(5,dom);
UU = diag(u(xx));
zero1 = norm( M - [ eye(15), zeros(15); DD, UU ])

%% more complicated chebmatrix
A = [ I, x, -3*I; 
    linop.sum(dom), 5, linop.feval(dom(end),dom);
    D, chebfun(1,dom), U ];
M = discretize(A,5);
MM = [ eye(15), xx, -3*eye(15);  
    ww, 5, [zeros(1,14) 1];
    DD, ones(15,1), UU ];
zero2 = norm( M - MM )

%%
% application to an appropriate chebmatrix
v = [ exp(x); pi; cos(x) ];
Av = A*v;
zero3 = norm( (v{1} + pi*x -3*v{3}) - Av{1} )
zero4 = norm( sum(v{1})+5*v{2}+feval(v{3},dom(end)) - Av{2} )
zero5 = norm( diff(v{1})+pi+(u).*v{3} - Av{3} )
