% Test for @chebfun/chebpade.m.
% Based on the Chebfun v4 test written by Ricardo Pachon, Feb. 09, 2009.

function pass = test_chebpade(pref)

% Turn off warnings and save the warning state.
warnState = warning('off');

% An example which doesn't require degree-reduction.
dom = [-1, 3];

P = chebfun([ 0.5045; -1.3813; 2.1122; 0.0558; -0.6817], dom, 'coeffs');
Q = chebfun([ 1; 0.1155; -0.8573; -0.2679; 0.5246 ], dom, 'coeffs');
R = P./Q;
[p, q, r] = chebpade(R, 4, 4, 'maehly');
err = norm(P - p) + norm(Q - q);
pass(1) = err < 100*vscale(R)*eps;

% An example which requires degree-reduction.
[p, q, r] = chebpade(P./Q, 6, 5, 'maehly');
err = norm(P - p) + norm(Q - q);
pass(2) = err < 100*vscale(R)*eps;

% An example by Geddes.
a = [-464/6375 ; -742/6375 ; 349/12750 ; 512/6375 ; 13/3400 ; 2129/51000 ; ...
     1333/8160 ; 9703/34000];
b = [-32/85 ; -28/85 ; 1];
f = chebfun(@(x) polyval(a, x)./polyval(b, x));
[p, q, r] = chebpade(f, 7, 2);
pass(3) = norm(f - p./q, inf) < 2e-15;

cp = chebcoeffs(p, length(p));
pass(4) = abs(cp(1) - 17/46) < 1e-13;

% Try non-rational and complex-valued functions.
f = chebfun(@exp);
[p, q, r] = chebpade(f, 2, 3);
[p2, q2, r2] = chebpade(1i*f, 2, 3);
ratio1 = p2(.3)/p(.3);
ratio2 = r2(-.4)/r(-.4);
pass(5) = norm([ratio1 ratio2] - [1i 1i], inf) < 1e-13;

% Restore the warning state.
warning(warnState);

end
