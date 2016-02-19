dom = [0,1];
x = chebfun('x',dom);
h = -1 + 0*x.^2;
N = chebop(@(x,u) h(x).*diff(u,2), dom);
N.lbc = @(u) u + 2*diff(u);
N.rbc = @(u) diff(u);

Nstar = adjoint(N);
