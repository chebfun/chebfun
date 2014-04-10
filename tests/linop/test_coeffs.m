function pass = test_coeffs
% TAD, 10 Jan 2014

%% Building blocks
dom = [ - 2  - 0.5 1 2];
I = operatorBlock.eye(dom);
D = operatorBlock.diff(dom);
Z = operatorBlock.zeros(dom);
x = chebfun('x', dom);
u = sin(x.^2);
U = operatorBlock.mult(u);   

%% Coefficient realizations
% Quasimatrix [ 1 ]
Ic = toCoeff(I);
err(1) = norm(Ic{1} - x.^0);
% Quasimatrix [ 1 0 ]
Dc = toCoeff(D);
err(2) = norm(Dc{1} - x.^0);
err(3) = norm(Dc{2} - 0);
% Quasimatrix [ x.^2 ]
Uc = toCoeff(U);
err(4) = norm(Uc{1} - u);

%%
A = D*(I+U*D);
Ac = toCoeff(A);
err(5) = norm(Ac{1} - u);
err(6) = norm(Ac{2} - (diff(u)+1));
err(7) = norm(Ac{3} - 0);

pass = err < 1e-15;

end
