function pass = test_operatorForm

%% Building blocks
dom = [-2 2];
I = linop.eye(dom);
D = linop.diff(dom);
Z = linop.zeros(dom);
x = chebfun('x', dom);
u = sin(x.^2);
U = linop.mult(u);   

%% Operator instantiations
eyeop = I.functionForm;
err(1) = norm( eyeop(u) - u );
diffop = D.functionForm;
TwoX = diffop(x.*x);   % should be the chebfun 2*x
TwoXAgain = D*(x.*x);
err(2) = norm( 2*x - TwoX );
multop = U.functionForm;
err(3) = norm( multop(x+1) - u.*(x+1) );   % should be zero
% shorthand notation, U*chebfun
err(4) = norm( multop(x+1) - U*(x+1) );

%%

pass = err < 1e-14;

end