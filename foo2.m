dom = [0,pi];
N = chebop(dom);
N.op = @(x,u) -diff(u,2);
N.lbc = @(u) [u;diff(u)];
N.rbc = @(u) [diff(u)];

[U,S,V] = svds(N,4);
diag(S)
subplot(2,1,1), plot(U)
subplot(2,1,2), plot(V)

norm(N*V - U*S)
