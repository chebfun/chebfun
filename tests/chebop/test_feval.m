function pass = test_feval(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

%% Scalar equation:

x = chebfun(@(x) x, [0 pi], pref);
u = sin(x);
Nu = diff(u, 2) + cos(u);
N = chebop(@(u) diff(u, 2) + cos(u));
err(1) = norm(feval(N, u) - Nu);
err(2) = norm(N(u) - Nu);

N = chebop(@(x, u) diff(u, 2) + cos(u));
err(3) = norm(feval(N, u) - Nu);
err(4) = norm(N(u) - Nu);
err(5) = norm(N*u - Nu);
err(6) = norm(feval(N, x, u) - Nu);
err(7) = norm(N(x, u) - Nu);

%% System of equations:

N = chebop(@(x, u, v) [diff(u, 2) + cos(v) ; diff(v, 2) - sin(u)]);
x = chebfun(@(x) x, [0 pi], pref);
u = sin(x);
v = exp(x);
uv = [u ; v];
Nuv = [diff(u, 2) + cos(v) ; diff(v, 2) - sin(u)];
err(8) = norm(feval(N, u, v) - Nuv);
err(9) = norm(feval(N, x, u, v) - Nuv);
err(10) = norm(feval(N, uv) - Nuv);
err(11) = norm(N(u, v) - Nuv);
err(12) = norm(N(x, u, v) - Nuv);
err(13) = norm(N(uv) - Nuv);
err(14) = norm(N*uv - Nuv);

%%

tol = 1e-14;
pass = err < tol;

end