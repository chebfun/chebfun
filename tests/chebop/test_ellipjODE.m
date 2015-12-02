function pass = test_ellipjODE(pref)
% Test to check that the chebfun command for Jacobi elliptic functions
% ELLIPJ produces a similar result to the solution of the nonlinear
% differential equation that generates it.

% NOTE: Taken from V4 test chebop_ellipjode.

if ( nargin == 0 )
    pref = cheboppref();
end
tol = 1e-10;

% Elliptic parameter
m = .1;
K = ellipke(m);
d = [0, K];
x = chebfun(@(x) x, d);

%% jacobi elliptic functions
[sn, cn, dn] = ellipj(x, m); 

% SN is the solution of an ODE
N = chebop(d);
N.op = @(x,u) diff(u,2) + (1+m)*u - 2*m*u.^3;
N.lbc = @(u) u;% The solution here is close to being singular, and we only enjoy linear convergence  so loosen the tolerance
N.rbc = @(u) u - 1 ;
u = N\0;
err(1) = norm(u - sn, inf);

%% CN
N.op = @(x, u) diff(u, 2) + (1-2*m)*u + 2*m*u.^3;
N.lbc = @(u) u - 1 ;
N.rbc = @(u) u ;
u = mldivide(N, 0, pref);
err(2) = norm(u - cn, inf);

% %% DN
% % The solution here is close to being singular, and we only enjoy linear
% % convergence of Newton iteratin. Loosen the tolerance.
% pref.bvpTol = 5e-5;
% N.op = @(x, u) diff(u, 2) - (2-m)*u + 2*u.^3;
% N.lbc = @(u) u - 1 ;
% N.rbc = @(u) u - sqrt(1-m);
% u = mldivide(N, 0, pref);
% err(3) = norm(u - dn, inf);

%%
% pass = err < [tol, tol, 2e6*tol];
pass = err < [tol, tol];

end


