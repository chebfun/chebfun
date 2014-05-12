function pass = test_eigenvaluesPiecewise(pref)
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
[Z, I, D, C, M] = linop.primitiveOperators(d);
[z, e, s, r] = linop.primitiveFunctionals(d);

% Build the operator
L = linop([ -h*D^2 + M(a*(sign(x-b)-sign(x-c)))]);
L = addbc( L, e(-5), 0 );
L = addbc( L, e(5), 0 );

% v4 results to compare against:
v4results = [0.055633547864; 
    0.058372224970;
    0.222492474974;
    0.233441822781;
    0.500447687842;
    0.525062886469];
%% Solve with colloc2
prefs = cheboppref;
prefs.discretization = @colloc2;
e = eigs(L, 6, prefs);
tol = 1e-10;
pass(1) = norm(e - v4results) < tol;

%% Solve with ultraS
prefs.discretization = @ultraS;
e = eigs(L, 6, 0, prefs);
tol = 1e-10;
pass(2) = norm(e - v4results, inf) < tol;

end
