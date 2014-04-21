function pass = test_linjump

% Test jump condition in an ODE.
% TAD, 31 Jan 2014

tol = 1e-9;

domain = [0 0.3 1];
[Z,I,D,C] = linop.primitiveOperators(domain);
[z,e,s] = linop.primitiveFunctionals(domain);
j = functionalBlock.jumpAt(domain);

A = linop( [ D^2, I; -D D^2+I ] );

A = addConstraint(A,[e(0) z],-1);
A = addConstraint(A,[z e(1)],1);
A = addContinuity(A,[j(0.3,1) z],2);
A = addContinuity(A,[j(0.3,0) j(0.3,0)],0);

x = chebfun('x',domain);
u = A\[x;0*x];

%%
% jumps
J = functionalBlock.jump(0.3,domain,0);
err(1) = J*u{1} + J*u{2};
err(2) = J*(D*u{1}) - 2;

%%
% BCs
err(3) = feval(u{1},0) + 1;
err(4) = feval(u{2},1) - 1;

%%
% ODEs
% [TODO]: The following lines use a hack. This needs to be changed:
g = D^2*u{1}+u{2};
g.funs{1} = g.funs{1}.funPart;
g.funs{2} = g.funs{2}.funPart;
err(5) = norm(g - x);

f = -D*u{1} + D^2*u{2} + u{2};
f.funs{1} = f.funs{1}.funPart;
f.funs{2} = f.funs{2}.funPart;
err(6) = norm(f);
%%
pass = abs(err) < tol;


end