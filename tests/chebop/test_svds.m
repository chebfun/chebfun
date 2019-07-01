function pass = test_svds(pref)

if ( nargin == 0 )
    pref = cheboppref();
end

tol = 1e1*pref.bvpTol;

%% SVD of differential operator
d = [0, pi];
N = chebop(d);
N.op = @(x,u) diff(u);
[U,S,V] = svds(N, 10);
s = diag(S);
s_true = (9:-1:0)';
pass(1) = norm(s-s_true,inf) < tol*max(s);
pass(2) = norm(N(V)-U*S) < tol*max(s);

end
