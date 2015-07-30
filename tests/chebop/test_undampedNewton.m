function pass = test_undampedNewton(pref)
%TEST_UNDAMPEDNEWTON    Test what happens when we turn off damping in Newton

% This test tries to solve steady state Allen-Cahn using undamped Newton
% iteration.

if ( nargin == 0 ) 
    pref = cheboppref();
end

% Turn off damping
pref.damping = 0;

%% Setup and solve
dom = [0, 10];
x = chebfun(@(x) x, dom);
f = sin(x);
N = chebop(@(u) diff(u,2) + u - u.^3, dom, 1, -1);
u = solvebvp(N, f, pref);

% Happy?
err = norm(N(u) - f);
pass = ( err < 1e-8 );

% end
