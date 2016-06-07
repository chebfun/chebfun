function pass = test_ellipke(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

%% 1 output:
% Linear scalar:
m = chebfun('m', [0 .99], pref);
K1 = ellipke(m);
K2 = chebfun(@(m) ellipke(m), [0 .99], pref);
pass(1) = normest(K1 - K2) < 1e2*eps*vscale(K1);


% Array-valued composition:
f = chebfun(@(x) .05+abs(.9*[sin(pi*x), cos(pi*x)]), -1:.5:1);
K1 = ellipke(f);
F = @(m) ellipke(.05 + abs(.9*[sin(pi*m), cos(pi*m)]));
K2 = chebfun(@(x) F(x), -1:.5:1, pref);
pass(2) = normest(K1 - K2) < 1e3*eps*vscale(K1);


%% 2 outputs:
% Linear scalar:
m = chebfun('m', [0 .99], pref);
[K1, E1] = ellipke(m);
E2 = chebfun(@(m) myellipke(m), [0 .99], pref);
err = normest(E1 - E2);
tol = 10*eps*vscale(E1);
pass(3) = err < tol;

% Array-valued composition:
f = chebfun(@(x) .05+abs(.9*[sin(pi*x), cos(pi*x)]), -1:.5:1);
[K1, E1] = ellipke(f);
F = @(m) myellipke(.05 + abs(.9*[sin(pi*m), cos(pi*m)]));
E2 = chebfun(@(x) F(x), -1:.5:1, pref);
pass(4) = normest(E1 - E2) < 1e3*eps*vscale(E1);


end

function E = myellipke(m)

[ignored, E] = ellipke(m);

end
