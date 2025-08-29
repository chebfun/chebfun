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

% Note. The existence of solutions to functional DEs is even more tricky
% that that of ODEs. Here we just demonstrate that the Chebfun 
% implementation at least seems to work in some cases. 

% Test 5 (periodic)
dom = [-pi, pi];
f = chebfun(@(x) exp(sin(x)), dom, 'periodic');
N = chebop(@(x,u) diff(u,2) - u(0) , dom, 'periodic');
u = N\f;
pass(5) = norm(N(u) - f) < tol;

% Test 6 (periodic)
N = chebop(@(x,u) diff(u,2) - u(x+2*pi) , dom, 'periodic');
u = N\f;
pass(6) = norm(N(u) - f) < tol;

% Test 7 (periodic)
N = chebop(@(x,u) diff(u,2) + u(u) , dom, 'periodic');
u = N\f;
pass(7) = norm(N(u) - f) < 10*tol;

end


