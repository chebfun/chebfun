function pass = test_eigenvaluesGeneralized(pref)
% Test generalized eigenvalueproblems for LINOP
if ( nargin == 0 )
    pref = cheboppref();
end

tol = 6e-9;

dom = [-1, 1];
D = operatorBlock.diff(dom);
ev = functionalBlock.eval(dom);

%%

% D^2*u = 1i*lam*D*u, u(-1) = u(1) = 0.
A = linop(D^2);
A = addbc(A, ev(-1));
A = addbc(A, ev(1));
B = linop(1i*D);

% True eigenvalues:
e_true = [-3 -2 -1 1 2 3].';

% COLLOC1
pref.discretization = @colloc1;
e = eigs(A, B, 6, 0, pref)/pi;
er = sort(real(e));
err(1,1) = norm(er - e_true) + norm(imag(e));

% COLLOC2
pref.discretization = @colloc2;
e = eigs(A, B, 6, 0, pref)/pi;
er = sort(real(e));
err(1,2) = norm(er - e_true) + norm(imag(e));

% ULTRAS
pref.discretization = @ultraS;
e = eigs(A, B, 6, 0, pref)/pi;
er = sort(real(e));
err(1,3) = norm(er - e_true) + norm(imag(e));


%%

% Now try putting the highest derivative on the right:

% 1i*D*u = (1/lam)*D*u, u(-1) = u(1) = 0.
A = linop(D^2);
B = linop(1i*D);
B = addbc(B, ev(-1));
B = addbc(B, ev(1));

% True eigenvalues:
e_true = (1:6).';

% COLLOC1
pref.discretization = @colloc1;
e = 1./eigs(B, A, 6, 1, pref)/pi;
er = sort(real(e));
err(2,1) = norm(er - e_true) + norm(imag(e));

% COLLOC2
pref.discretization = @colloc2;
e = 1./eigs(B, A, 6, 1, pref)/pi;
er = sort(real(e));
err(2,2) = norm(er - e_true) + norm(imag(e));

% ULTRAS
pref.discretization = @ultraS;
e = 1./eigs(B, A, 6, 1, pref)/pi;
er = sort(real(e));
err(2,3) = norm(er - e_true) + norm(imag(e));

%%

% Now try a sytem:

dom = [-1, 1];
D = operatorBlock.diff(dom);
I = operatorBlock.eye(dom);
Z = operatorBlock.zeros(dom);
ev = functionalBlock.eval(dom);
z = functionalBlock.zero(dom);

A = [D^2 D ; D D^2];
B = [I D ; D I];
A = linop(A);
A = addbc(A, [ev(-1) z]);
A = addbc(A, [ev(1) z]);
A = addbc(A, [z ev(-1)]);
A = addbc(A, [z ev(1)]);
B = linop(B);
e_true = -1 + 1i*pi*[-1 1 -1 1 -2 2 -2 2 -3 3 -3 3].';

% COLLOC1
pref.discretization = @colloc1;
e = eigs(A, B, 12, 0, pref);
err(3,1) = norm(e - e_true);

% COLLOC2
pref.discretization = @colloc2;
e = eigs(A, B, 12, 0, pref);
err(3,2) = norm(e - e_true);

% ULTRAS
pref.discretization = @ultraS;
e = eigs(A, B, 12, 0, pref);
err(3,3) = norm(e - e_true);

%%

pass = err < tol;

end