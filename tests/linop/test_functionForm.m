function pass = test_functionForm

%% Building blocks
dom = [-2 2];
I = operatorBlock.eye(dom);
D = operatorBlock.diff(dom);
Z = operatorBlock.zeros(dom);
x = chebfun('x', dom);
u = sin(x.^2);
U = operatorBlock.mult(u);   

%% Operator instantiations
eyeop = toFunction(I);
err(1) = norm( eyeop(u) - u );
diffop = toFunction(D);
TwoX = diffop(x.*x);   % should be the chebfun 2*x
TwoXAgain = D*(x.*x);
err(2) = norm( 2*x - TwoX );
multop = toFunction(U);
err(3) = norm( multop(x+1) - u.*(x+1) );   % should be zero
% shorthand notation, U*chebfun
err(4) = norm( multop(x+1) - U*(x+1) );

%%

pass = err < 1e-14;

end
