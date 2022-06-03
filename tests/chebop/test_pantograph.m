function pass = test_pantograph(pref)

tol = 1e-10;

if ( nargin == 0 )
    pref = cheboppref();
end
pref.discretization = 'values';

% Test 1:
x = chebfun('x', [0 2]);
N = chebop(@(x,u) diff(u,2) - u(x/2), [0 2]);
N.lbc = 0; N.rbc = 1;
u = N\0;
pass(1) = norm(diff(u,2) - u(x/2)) < tol;

% Test 2:
N = chebop(@(x,u) diff(u,2) - u(u), [0 2]);
N.lbc = 0; N.rbc = 1;
u = N\0;
pass(2) = norm(diff(u,2) - u(u)) < tol;

% Test 3:
x = chebfun('x', [0 2]);
N = chebop(@(x,u) diff(u,2) - u(.5*diff(u)), [0 2]);
N.lbc = 0; N.rbc = 1;
u = N\0;
pass(3) = norm(diff(u,2) - u(.5*diff(u))) < tol;

% Test 4
N = chebop(@(x,u) diff(u,2) - u(x*u^2) , [0 2]);
N.lbc = 0.1; N.rbc = 1;
u = N\0;
pass(4) = norm(diff(u,2) - u(x*u^2)) < 100*tol;