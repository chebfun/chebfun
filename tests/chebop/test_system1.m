function pass = test_system1(pref)
% Test 2x2 system (sin/cos)
% Toby Driscoll

% Note: This test is taken from chebop_systemsolve1 in V4.

if ( nargin == 0 )
    pref = cheboppref;
end
pref.plotting = 'off';
pref.display = 'off';

tol = 1e-10;

% Smooth domain:
d = [-pi pi];
A = chebop(@(x,u,v) [u - diff(v) ; diff(u) + v],d);
A.lbc = @(u,v) u + 1;
A.rbc = @(u,v) v;
x = chebfun('x',d);
f = [ 0*x ; 0*x ];
u = mldivide(A, f, pref);

u1 = u{1}; u2 = u{2};
pass(1) = norm( u1 - cos(x),inf) < 100*tol;
pass(2) = norm( u2 - sin(x),inf) < 100*tol;

%% Piecewise:
A.domain = [-pi 0 pi];
u = mldivide(A, f, pref);
u1 = u{1}; u2 = u{2};
pass(3) = norm( u1 - cos(x),inf) < 2000*tol;
pass(4) = norm( u2 - sin(x),inf) < 2000*tol;

end


