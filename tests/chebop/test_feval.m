function pass = test_feval(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

%% Scalar equation:
dom = [0 pi];
x = chebfun(@(x) x, dom, pref);
u = sin(x);
Nu = diff(u, 2) + cos(u);
N = chebop(@(u) diff(u, 2) + cos(u), dom);
err(1) = norm(feval(N, u) - Nu);
err(length(err) + 1) = norm(N(u) - Nu);

N = chebop(@(x, u) diff(u, 2) + cos(u), dom);
err(length(err) + 1) = norm(feval(N, u) - Nu);
err(length(err) + 1) = norm(N(u) - Nu);
err(length(err) + 1) = norm(N*u - Nu);
err(length(err) + 1) = norm(feval(N, x, u) - Nu);
err(length(err) + 1) = norm(N(x, u) - Nu);

%% Scalar equation, chebmatrix input:
N = chebop(@(u) diff(u, 2) + cos(u), dom);
u = chebmatrix({u});

err(length(err) + 1) = norm(feval(N, u) - Nu);
err(length(err) + 1) = norm(N(u) - Nu);

N = chebop(@(x, u) diff(u, 2) + cos(u), dom);
err(length(err) + 1) = norm(feval(N, u) - Nu);
err(length(err) + 1) = norm(N(u) - Nu);
err(length(err) + 1) = norm(N*u - Nu);
err(length(err) + 1) = norm(feval(N, x, u) - Nu);
err(length(err) + 1) = norm(N(x, u) - Nu);
%% System of equations:

N = chebop(@(x, u, v) [diff(u, 2) + cos(v) ; diff(v, 2) - sin(u)], dom);
x = chebfun(@(x) x, dom, pref);
u = sin(x);
v = exp(x);
uv = [u ; v];
Nuv = [diff(u, 2) + cos(v) ; diff(v, 2) - sin(u)];
err(length(err) + 1) = norm(feval(N, u, v) - Nuv);
err(length(err) + 1) = norm(feval(N, x, u, v) - Nuv);
err(length(err) + 1) = norm(feval(N, uv) - Nuv);
err(length(err) + 1) = norm(N(u, v) - Nuv);
err(length(err) + 1) = norm(N(x, u, v) - Nuv);
err(length(err) + 1) = norm(N(uv) - Nuv);
err(length(err) + 1) = norm(N*uv - Nuv);

%% Quasimatrix notation, part I:
N = chebop(@(x, u) [diff(u{1}, 2) + cos(u{2}) ; diff(u{2}, 2) - sin(u{1})]);
x = chebfun(@(x) x, [0 pi], pref);
u = [sin(x) ; exp(x)];
Nu = [diff(u{1}, 2) + cos(u{2}) ; diff(u{2}, 2) - sin(u{1})];
err(length(err) + 1) = norm(feval(N, u) - Nu);
err(length(err) + 1) = norm(feval(N, x, u) - Nu);
err(length(err) + 1) = norm(N(u) - Nu);
err(length(err) + 1) = norm(N(x, u) - Nu);
err(length(err) + 1) = norm(N*u - Nu);

%% Quasimatrix notation, part II:
N = chebop(@(u) [diff(u{1}, 2) + cos(u{2}) ; diff(u{2}, 2) - sin(u{1})]);
x = chebfun(@(x) x, [0 pi], pref);
u = [sin(x) ; exp(x)];
Nu = [diff(u{1}, 2) + cos(u{2}) ; diff(u{2}, 2) - sin(u{1})];
err(length(err) + 1) = norm(feval(N, u) - Nu);
err(length(err) + 1) = norm(N(u) - Nu);
err(length(err) + 1) = norm(N*u - Nu);

%% Test multiple RHS on an e-val problem: (From #1022)
dom = [0 pi];
op = @(x,u,v) [ -diff(v) ; diff(u) ];
bc = @(x,u,v) [ u(0) ; u(pi) ];
A = chebop(op, dom, bc);
[ev, ew] = eigs(A);
err(length(err) + 1) = norm(A*ev - ev*ew)/1e3;

%% Happy?

tol = 1e-14;
pass = err < tol;

end
