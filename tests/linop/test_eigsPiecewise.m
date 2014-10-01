function pass = test_eigsPiecewise(pref)
% TEST_EIGVENVALUESPIECEWISE    Solve Schrodinger equation with a piecewise well
% AB + H, 2014/03/10

%% Setup

if ( nargin == 0 )
    pref = chebfunpref();
end

% Domain, RHS, constants
d = [-5,5];
x = chebfun('x',d);
h = 0.1;
a = 4;
b = -1;
c = 0.9;

% Primitive operator blocks
[Z, I, diffOp, C, M] = linop.primitiveOperators(d);
[z, e, s, r] = linop.primitiveFunctionals(d);

% Build the operator
L = linop(-h*diffOp^2 + M(a*(sign(x-b)-sign(x-c))));
L = addbc( L, e(-5), 0 );
L = addbc( L, e(5), 0 );

% v4 results to compare against:
v4results = [0.055633547864; 
    0.058372224970;
    0.222492474974;
    0.233441822781;
    0.500447687842;
    0.525062886469];

prefs = cheboppref;

%% Solve with chebcolloc1
prefs.discretization = @chebcolloc1;
[V, D] = eigs(L, 6, prefs);
e = diag(D);
err(1) = norm(e - v4results, inf);
% Check that we actually computed eigenfunctions
err(2) = norm(L*V-V*D);
%% Solve with chebcolloc2
prefs.discretization = @chebcolloc2;
[V, D] = eigs(L, 6, prefs);
e = diag(D);
err(3) = norm(e - v4results, inf);
% Check that we actually computed eigenfunctions
err(4) = norm(L*V-V*D);

%% Solve with ultraS
prefs.discretization = @ultraS;
[V, D] = eigs(L, 6, prefs);
e = diag(D);
err(5) = norm(e - v4results, inf);
% Check that we actually computed eigenfunctions
err(6) = norm(L*V-V*D);

%%
tol = [1e-10 5e-8 1e-10 5e-8 1e-10 5e-8];
pass = err < tol;

end
