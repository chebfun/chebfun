function pass = test_mtimes(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

% Note that N*u where u is a chebfun or chebmatrix is tested in test_feval.

%% Scalar equation:

x = chebfun(@(x) x, [0 pi], pref);
u = sin(x);
a = sqrt(2);
N = chebop(@(u) diff(u, 2) + cos(u));
err(1) = norm((a*N)*u - a*(N*u));
err(2) = norm((N*a)*u - a*(N*u));


%% System of equations:

N = chebop(@(x, u, v) [diff(u, 2) + cos(v) ; diff(v, 2) - sin(u)]);
x = chebfun(@(x) x, [0 pi], pref);
u = sin(x);
v = exp(x);
uv = [u ; v];
a = sqrt(2);
err(3) = norm((a*N)*uv - a*(N*uv));
err(4) = norm((N*a)*uv - a*(N*uv));

%% Some things that shouldn't work:

try
    aaN = [a a]*N;
    err(5) = inf;
catch
    err(5) = 0;
end

try
    NN = N*N;
    err(6) = inf;
catch
    err(6) = 0;
end

try
    uN = u*N;
    err(7) = inf;
catch
    err(7) = 0;
end

%%

tol = 1e-14;
pass = err < tol;

end
