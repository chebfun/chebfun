function pass = test_eigenvaluesGeneralized(pref)
% Test generalized eigenvalueproblems for LINOP
if ( nargin == 0 )
    pref = cheboppref();
end

dom = [-1, 1];
diffOp = operatorBlock.diff(dom);
ev = functionalBlock.eval(dom);

%%

% D^2*u = 1i*lam*D*u, u(-1) = u(1) = 0.
A = linop(diffOp^2);
A = addbc(A, ev(-1));
A = addbc(A, ev(1));
B = linop(1i*diffOp);

% True eigenvalues:
e_true = [-3 -2 -1 1 2 3].';

% COLLOC1
pref.discretization = @colloc1;
[V, D] = eigs(A, B, 6, 0, pref);
e = diag(D)/pi;
er = sort(real(e));
AV = A*V;
BV = B*V;
err(1,1) = norm(er - e_true) + norm(imag(e));
err(1,2) = norm(AV-BV*D);

% COLLOC2
pref.discretization = @colloc2;
[V, D] = eigs(A, B, 6, 0, pref);
e = diag(D)/pi;
er = sort(real(e));
AV = A*V;
BV = B*V;
err(1,3) = norm(er - e_true) + norm(imag(e));
err(1,4) = norm(AV-BV*D);

% ULTRAS
pref.discretization = @ultraS;
[V, D] = eigs(A, B, 6, 0, pref);
e = diag(D)/pi;
er = sort(real(e));
AV = A*V;
BV = B*V;
err(1,5) = norm(er - e_true) + norm(imag(e));
err(1,6) = norm(AV-BV*D);



%%

% Now try putting the highest derivative on the right:

% 1i*D*u = (1/lam)*D*u, u(-1) = u(1) = 0.
A = linop(diffOp^2);
B = linop(1i*diffOp);
B = addbc(B, ev(-1));
B = addbc(B, ev(1));

% True eigenvalues:
e_true = (1:6).';

% COLLOC1
pref.discretization = @colloc1;
[V, D] = eigs(B, A, 6, 1, pref);
e = 1./diag(D)/pi;
er = sort(real(e));
AV = A*V;
BV = B*V;
err(2,1) = norm(er - e_true) + norm(imag(e));
err(2,2) = norm(BV-AV*D);

% COLLOC2
pref.discretization = @colloc2;
[V, D] = eigs(B, A, 6, 1, pref);
e = 1./diag(D)/pi;
er = sort(real(e));
AV = A*V;
BV = B*V;
err(2,3) = norm(er - e_true) + norm(imag(e));
err(2,4) = norm(BV-AV*D);

% ULTRAS
pref.discretization = @ultraS;
[V, D] = eigs(B, A, 6, 1, pref);
e = 1./diag(D)/pi;
er = sort(real(e));
AV = A*V;
BV = B*V;
err(2,5) = norm(er - e_true) + norm(imag(e));
err(2,6) = norm(BV-AV*D);

%%

% Now try a sytem:

dom = [-1, 1];
diffOp = operatorBlock.diff(dom);
I = operatorBlock.eye(dom);
Z = operatorBlock.zeros(dom);
ev = functionalBlock.eval(dom);
z = functionalBlock.zero(dom);

A = [diffOp^2 diffOp ; diffOp diffOp^2];
B = [I diffOp ; diffOp I];
A = linop(A);
A = addbc(A, [ev(-1) z]);
A = addbc(A, [ev(1) z]);
A = addbc(A, [z ev(-1)]);
A = addbc(A, [z ev(1)]);
B = linop(B);
e_true = -1 + 1i*pi*[-1 1 -1 1 -2 2 -2 2 -3 3 -3 3].';

% COLLOC1
pref.discretization = @colloc1;
[V, D] = eigs(A, B, 12, 0, pref);
e = diag(D);
AV = A*V;
BV = B*V;
err(3,1) = norm(e - e_true);
err(3,2) = norm(AV - BV*D);

% COLLOC2
pref.discretization = @colloc2;
[V, D] = eigs(A, B, 12, 0, pref);
e = diag(D);
AV = A*V;
BV = B*V;
err(3,3) = norm(e - e_true);
err(3,4) = norm(AV - BV*D);

% ULTRAS
pref.discretization = @ultraS;
[V, D] = eigs(A, B, 12, 0, pref);
e = diag(D);
AV = A*V;
BV = B*V;
err(3,5) = norm(e - e_true);
err(3,6) = norm(AV - BV*D);

%%
tolVals = repmat(6e-9, 3, 1);
tolFuns = repmat(4e-7, 3, 1);

tol = repmat([tolVals, tolFuns], 1, 3);

pass = err < tol;

end
