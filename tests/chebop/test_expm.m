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

err(1) = norm(V4 - feval(u{6}, chebpts(15)));

%%
% Test backward compatability:
A = chebop(@(u) diff(u, 2), [-1, 1], 0);
E = expm(A, t);
u6 = E*u0;

err(2) = norm(V4 - feval(u6, chebpts(15)));

%%

pass = err < tol;

end