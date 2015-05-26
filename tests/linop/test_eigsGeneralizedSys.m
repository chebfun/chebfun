function pass = test_eigsGeneralizedSys(pref)
% Test generalized eigenvalueproblems for LINOP
if ( nargin == 0 )
    pref = cheboppref();
end

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

% CHEBCOLLOC1
pref.discretization = @chebcolloc1;
[V, D] = eigs(A, B, 12, 0, pref);
e = diag(D);
AV = A*V;
BV = B*V;
err(1,1) = norm(e - e_true);
err(1,2) = norm(AV - BV*D);

% CHEBCOLLOC2
pref.discretization = @chebcolloc2;
[V, D] = eigs(A, B, 12, 0, pref);
e = diag(D);
AV = A*V;
BV = B*V;
err(1,3) = norm(e - e_true);
err(1,4) = norm(AV - BV*D);

% ULTRAS
pref.discretization = @ultraS;
[V, D] = eigs(A, B, 12, 0, pref);
e = diag(D);
AV = A*V;
BV = B*V;
err(1,5) = norm(e - e_true);
err(1,6) = norm(AV - BV*D);

%%
tolVals = 6e-9;
tolFuns = 4e-7;

tol = repmat([tolVals, tolFuns], 1, 3);

pass = err < tol;

end
