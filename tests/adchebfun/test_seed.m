function pass = test_seed(pref)
% Test for ADCHEBFUN/SEED()
% AB, 2014/05/14

if ( nargin == 0 )
    pref = chebfunpref();
end

%% Setup
dom = [-2 4];
x = chebfun(@(x) x, dom);
u = adchebfun(x);
tol = 1e-15;
%% Try various calls to seed
% Seed as identity operator
u = seed(u, 1, 1);
uJac = u.jacobian;
err =  norm(feval(toFunction(uJac), x) - x);
pass(1) = ( err < tol );

% Seed as identity function
u = seed(u, 1, 0);
uJac = u.jacobian;
err = norm(uJac - 1);
pass(2) = ( err < tol );

% Shortcut seeding
u = seed(u, 2, 3);
uJac = u.jacobian;
err = norm(feval(toFunction(uJac{2}), x) - x) + ...
    norm(feval(toFunction(uJac{1}), x)) + norm(feval(toFunction(uJac{3}), x));
pass(3) = ( err < tol );

% Combination of operators and functions
boolVec = [true, false, false, true];
u = seed(u, 3, boolVec);
uJac = u.jacobian;
err = norm(feval(toFunction(uJac{1}), x)) + norm(uJac{2}) + ...
    norm(uJac{3} - 1) + norm(feval(toFunction(uJac{4}), x));
pass(4) = ( err < tol );

% Swap roles!
u = seed(u, 3, ~boolVec);
uJac = u.jacobian;
err = norm(uJac{1}) + norm(feval(toFunction(uJac{2}), x)) + ...
    norm(feval(toFunction(uJac{3}), x) - x) + norm(uJac{4});
pass(5) = ( err < tol );
