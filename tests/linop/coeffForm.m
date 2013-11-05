clc
clear classes

%% Building blocks
dom = [-2 -0.5 1 2];
I = linop.eye(dom);
D = linop.diff(dom);
Z = linop.zeros(dom);
x = chebfun('x', dom);
u = sin(x.^2);
U = linop.mult(u);   

%% Coefficient realizations
% Quasimatrix [ 1 ]
Ic = I.coeffForm;
zero1 = norm( Ic{1} - x.^0 )
% Quasimatrix [ 1 0 ]
Dc = D.coeffForm;
zero2 = norm( Dc{1} - x.^0 )
zero3 = norm( Dc{2} - 0 )
% Quasimatrix [ x.^2 ]
Uc = U.coeffForm;
zero4 = norm( Uc{1} - u )

%%
A = D*(I+U*D);
Ac = A.coeffForm;
zero5 = norm( Ac{1}-u )
zero6 = norm( Ac{2}-(diff(u)+1) )
zero7 = norm( Ac{3}-0 )