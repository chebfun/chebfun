rng(1);
dom = [-pi,pi];
x = chebfun('x',dom);
a3 = 1+ sin(x).^2; a2 = 1 + x.^2; a1 = sin(x); a0 = exp(x);
br2 = randn; br1 = randn; br0 = randn;
bl12 = randn; bl11 = randn; bl10 = randn;
bl22 = randn; bl21 = randn; bl20 = randn;
N = chebop(@(x,u) a3.*diff(u,3) + a2.*diff(u,2) + a1.*diff(u) + a0.*u, dom);
lbcOp = @(u) [ bl12*diff(u,2) + bl11*diff(u) + bl10*u;...
               bl22*diff(u,2) + bl21*diff(u) + bl20*u ];
rbcOp = @(u) br2*diff(u,2) + br1*diff(u) + br0*u;
N.lbc = lbcOp; 
N.rbc = rbcOp;

Nstar = adjoint(N)

f = chebpoly(5,dom);
u = N\f; 
v = Nstar\f;
plot([f,u,v])

i1 = v'*(N*u);
i2 = (Nstar*v)'*u;

abs(i1-i2)/abs(i1)
