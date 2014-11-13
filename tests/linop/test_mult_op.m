function pass = test_mult_op
% Check if multiplication with a pointwise multiplication operator
% is identical to a pointwise multiplication of two chebfuns.

% TAD, 3 Feb 2014

tol = 1e-14;

d = [0,2];
x = chebfun('x',d);
f = sin(exp(2*x));
g = x.^3-cos(x);

% no breaks
F = operatorBlock.mult(f);
err(1) = norm( F*g - f.*g );

% breakpoint
f = abs(f);
F = operatorBlock.mult(f);
err(2) = norm( F*g - f.*g );

%% Expansion
% Try to let an multiplication operator operate on an array-valued CHEBFUN:
h = cos(x);
M = operatorBlock.mult(h);
V = chebpoly(1:4, d);
err(3) = norm( M*V - repmat(h, 1, 4).*V);

%%
pass = err < tol;


end
