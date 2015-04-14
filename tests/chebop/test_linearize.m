function pass = test_linearize(pref)

if ( nargin == 0 )
    pref = cheboppref();
end
tol = 1e-14;

% TODO: More extensive testing is required.

%% Some toy problems:

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

% Pass a flag that we don't want to detect parameters:
N = chebop(@(u) u.^2, dom);
L = linearize(N, u, [], false, false);
err(3) = norm(L*v - 2*u.*v);

N = chebop(@(u) diff(u,2) + u.^2, dom);
L = linearize(N, u);
err(4) = norm(L*v - (diff(v, 2) + 2*u.*v));

N = chebop(@(x, u, v) [diff(u,2) + v.^2 ; diff(v) - cos(u)], dom);
L = linearize(N, [u ; v]);
err(5) = norm(L*[w ; v] - [diff(w, 2) + 2*v.^2 ; diff(v) + sin(u).*w]);

%% Test some problems with parameters:

N = chebop(@(x, u, a) diff(u,2) + a.^2, dom);
u0 = [chebfun(0, dom) ; 1];
L = linearize(N, u0);   % With init
err(6) = ~isa(L.blocks{2}, 'chebfun');
L = linearize(N);       % Without init
err(7) = ~isa(L.blocks{2}, 'chebfun');

%% Coupled system with four unknown variables

% Define the chebop
N = chebop(@(x,u,v,a,b) [x.*v + .001*diff(u,2) + a + 2*b ; diff(a.*v)-u], [-1 1]);
% First, we explicitly pass a chebmatrix to linearize around that contains two
% chebfuns and two doubles. The code should respect those dimensions, even
% though we encounter diff(a.*v) (which otherwise indicate that a should be
% treated as a function):
u0 = [chebfun(@sin) ; chebfun(@cos) ; 1 ; 0];
L = linearize(N, u0);   % With init
tmp = cellfun(@(b) isa(b, 'chebfun'), L.blocks) ~= [0 0 1 1 ; 0 0 1 1];
err(8) = any(tmp(:));

% If we don't pass a chebmatrix to linearize around, the code should assume that
% a is to be treated as a function (as it appears in diff(a.*v)). Hence, the
% dimensions of the output derivative will be different from before.
L = linearize(N);       % Without init
tmp = cellfun(@(b) isa(b, 'chebfun'), L.blocks) ~= [0 0 0 1 ; 0 0 0 1];
err(9) = any(tmp(:));

%%

pass = err < tol;

end
