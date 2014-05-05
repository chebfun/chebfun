function pass = test_linearize(pref)

if ( nargin == 0 )
    pref = cheboppref();
end
tol = 1e-14;

% TODO: More extensive testing is required.
% TODO: Support for parameterized problems is required.

dom = [0, 1, pi];
x = chebfun('x', dom);
u = sin(x);
v = cos(x);
w = exp(-x);

N = chebop(@(u) diff(u), dom);
L = linearize(N);
err(1) = norm(L*u - diff(u));

N = chebop(@(u) diff(u, 2) + diff(u) + sin(x).*u, dom);
L = linearize(N);
err(2) = norm(L*u - (diff(u, 2) + diff(u) + sin(x).*u));

N = chebop(@(u) u.^2, dom);
L = linearize(N, u);
err(3) = norm(L*v - 2*u.*v);

N = chebop(@(u) diff(u,2) + u.^2, dom);
L = linearize(N, u);
err(4) = norm(L*v - (diff(v, 2) + 2*u.*v));

N = chebop(@(x, u, v) [diff(u,2) + v.^2 ; diff(v) - cos(u)], dom);
L = linearize(N, [u ; v]);
err(5) = norm(L*[w ; v] - [diff(w, 2) + 2*v.^2 ; diff(v) + sin(u).*w]);

pass = err < tol;

end
