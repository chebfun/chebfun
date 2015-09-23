function pass = test_expm(pref)

if ( nargin == 0 )
    pref = cheboppref();
end

tol = pref.errTol;

%%
% Test against a V4 computation:

A = chebop(@(u) diff(u, 2), [-1, 1], 0);
x = chebfun('x');
u0 = exp(-20*(x+0.3).^2);  
t = [0 0.001 0.01 0.1 0.5 1];
u = expm(A, t, u0, pref);

V4 = [   ...
	0.000000000000000
    0.001144467788672
    0.004503386286796
    0.009766928773117
    0.016203479337652
    0.022565344171707
    0.027291616930494
    0.029038055428660
    0.027272749752324
    0.022536655030494
    0.016176338870674
    0.009748370213284
    0.004494411027930
    0.001142160800217
    0.000000000000000];

pass(1) = norm(V4 - feval(u{6}, chebpts(15))) < tol;

%%
% Test backward compatibility:
A = chebop(@(u) diff(u, 2), [-1, 1], 0);
warnState = warning('off', 'CHEBFUN:CHEBOP:expm:deprecated');
E = expm(A, t);
warning(warnState);
u6 = E*u0;

pass(2) = norm(V4 - feval(u6, chebpts(15))) < tol;

%%
% Test periodic boundary conditions with TRIGCOLLOC.
dom = [0 2*pi];
A = chebop(@(u) diff(u, 2), dom);
A.bc = 'periodic';
u0 = chebfun(@(x) sin(x), dom); 
t = [0 0.001 0.01 0.1 0.5 1];

% Solve with FOURIER technology in value space.
u = expm(A, t, u0);
pass(3) = isequal(get(u{1}.funs{1}, 'tech'), @trigtech);

% Solve with CHEBYSHEV technology in value space.
pref.discretization = @chebcolloc2;
v = expm(A, t, u0, pref);
pass(4) = isequal(get(v{1}.funs{1}, 'tech'), @chebtech2);

% Solve with CHEBYSHEV technology in coefficient space.
pref.discretization = @ultraS;
w = expm(A, t, u0, pref);
pass(5) = isequal(get(w{1}.funs{1}, 'tech'), @chebtech2);

% Compare solutions at final time.
pass(6) = norm(u{6} - v{6}, inf) < tol;
pass(7) = norm(w{6} - v{6}, inf) < tol;



end
