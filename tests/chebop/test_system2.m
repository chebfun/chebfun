function pass = test_system2(pref)
% Test solution of a 2x2 system
% Toby Driscoll

% Note: This test is taken from chebop_systemsolve1 in V4.

if ( nargin == 0 )
    pref = cheboppref;
end
pref.plotting = 'off';
pref.display = 'off';

tol = 1e-10;

% Smooth domain:
d = [-1 1];
A = chebop(@(x,u,v) [diff(u) + u + 2*v ; diff(u) - u + diff(v)], d);
A.lbc = @(u,v) u+diff(u);
A.rbc = @(u,v) diff(v);
x = chebfun('x',d);
f = [ exp(x) ; chebfun(1,d) ];
u = mldivide(A, f, pref);

u1 = u{1}; u2 = u{2};
pass(1) = norm( diff(u1)+u1+2*u2-exp(x)) < tol;
pass(2) = norm( diff(u1)-u1+diff(u2)-1 ) < tol;

%% Piecewise:
A.domain = [-1 0 1];
u = mldivide(A, f, pref);
u1 = u{1}; u2 = u{2};

err1 = diff(u1) + u1 + 2*u2 - exp(x);
err2 = diff(u1) - u1 + diff(u2) - 1;

pass(3) = norm( err1 ) < tol;
pass(4) = norm( err2 ) < tol;

end
