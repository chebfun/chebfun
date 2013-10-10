clc
clear classes

%% Building blocks
dom = [-2 2];
I = linop.eye(dom);
D = linop.diff(dom);
Z = linop.zeros(dom);
x = chebfun('x', dom);
u = sin(x.^2);
U = linop.diag(u);   

%% Operator instantiations
eyeop = op(I);
zero1 = norm( eyeop(u) - u )
diffop = op(D);
TwoX = diffop(x.*x);   % should be the chebfun 2*x
TwoXAgain = D*(x.*x);
zero2 = norm( 2*x - TwoX )
multop = op(U);
zero3 = norm( multop(x+1) - u.*(x+1) )   % should be zero
% shorthand notation, U*chebfun
zero4 = norm( multop(x+1) - U*(x+1) )
