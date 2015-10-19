function pass = test_diag(pref)
% Test for CHEBFUN/DIAG()

if ( nargin < 1 )
    pref = chebfunpref();
end

%% Test a few things
% No breakpoints
x = chebfun(@(x) x, [-1 6]);
tol = 1e-15;
f = sin(x);
g = cos(4*x);

D = diag(f);
err = norm(D*g -f.*g);
pass(1) = ( err < tol );

%% Breakpoints in g
x = chebfun(@(x) x, [-1 4 6]);
g = cos(4*x);
err = norm(D*g -f.*g);
pass(2) = ( err < tol );

%% Breakpoints in f and g
x = chebfun(@(x) x, [-1 4 6]);
f = sin(x);
g = cos(4*x);
err = norm(D*g -f.*g);
pass(3) = ( err < 10*tol );

end
