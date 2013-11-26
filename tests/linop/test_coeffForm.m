function pass = test_coeffForm


%% Building blocks
dom = [ - 2  - 0.5 1 2];
I = linop.eye(dom);
D = linop.diff(dom);
Z = linop.zeros(dom);
x = chebfun('x', dom);
u = sin(x.^2);
U = linop.mult(u);   

%% Coefficient realizations
% Quasimatrix [ 1 ]
Ic = I.coeffForm;
err(1) = norm(Ic{1} - x.^0);
% Quasimatrix [ 1 0 ]
Dc = D.coeffForm;
err(2) = norm(Dc{1} - x.^0);
err(3) = norm(Dc{2} - 0);
% Quasimatrix [ x.^2 ]
Uc = U.coeffForm;
err(4) = norm(Uc{1} - u);

%%
A = D*(I+U*D);
Ac = A.coeffForm;
err(5) = norm(Ac{1} - u);
err(6) = norm(Ac{2} - (diff(u)+1));
err(7) = norm(Ac{3} - 0);

err;
pass = err < 1e-15;

end
