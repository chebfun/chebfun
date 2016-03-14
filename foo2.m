dom = [0,pi];
N = chebop(dom);
N.op = @(x,u) -diff(u,2);
N.bc = 'dirichlet';


[U,S,V] = svds(N,4,0);
diag(S)
subplot(2,1,1), plot(U)
subplot(2,1,2), plot(V)
