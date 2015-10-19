% Test file for fov.m.
% Based on the test from Chebfun v4 written by Rodrigo Platte, Feb. 2009.

function pass = test_fov(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

A = [1 2 ; 3 2i];
F = fov(A);
err = abs(max(real(F)) - 3.049509756796393);
tol = 1e2*vscale(F)*eps;
pass(1) = err < tol;


B = 7;
F = fov(B);
err = norm(abs(F - 7));
tol = 10*vscale(F)*eps;
pass(2) = err < tol;

C = diag([-1 1 1i]);
F = fov(C);
err = abs(mean(F) - .25i);
tol = 10*vscale(F)*eps;
pass(3) = err < tol;

A = [0 1 0 ; 0 0 0 ; 0 0 1];
[F, lineSegs, theta] = fov(A);
err = norm(theta - pi/3*[1 5], inf);
pass(4) = err < tol;
tmp = feval(lineSegs, [-1;1]);
err = norm(real(tmp) - [1 0.25 ; 0.25 1], inf);
pass(5) = err < tol;

end
