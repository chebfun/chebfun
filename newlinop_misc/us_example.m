ccc

N = .001*linop.diff([-1, 1], 2) + linop.diff + linop.eye([-1, 1]); 
B = chebmatrix({linop.feval(-1, [-1 1])});
B2 = chebmatrix({linop.feval(1, [-1 1])});

A = chebmatrix({N});
L = linop(A); 
bc = linopConstraint();
bc = append(bc, B, 1);
bc = append(bc, B2, 1);
L.constraint = bc;

% discretize(L.operator, 5)
% discretize(L.constraint.operator, 5)

f = chebfun(0);

u = linsolve(L,f); 
u = u{1}
plot(u), shg