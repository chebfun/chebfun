% Test file for fov.m.
% Based on the test from Chebfun v4 written by Rodrigo Platte, Feb. 2009.

function pass = test_fov(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

A = [1 2 ; 3 2i];
F = fov(A);
err = abs(max(real(F)) - 3.049509756796393);
tol = 10*vscale(F)*epslevel(F);
pass(1) = err < tol;

B = 7;
F = fov(B);
err = norm(abs(F - 7));
tol = 10*vscale(F)*epslevel(F);
pass(2) = err < tol;

C = diag([-1 1 1i]);
F = fov(C);
err = abs(mean(F) - .25i);
tol = 10*vscale(F)*epslevel(F);
pass(3) = err < tol;

end
