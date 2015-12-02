function pass = test_paramODE_nonlin_C1(pref)
% Test solving a parameter dependent ODE. 
% Nick Hale, August 2011

% NOTE: Taken from V4 test chebop_paramODE.

if ( nargin == 0 )
    pref = cheboppref();
end
tol = 1e3*pref.bvpTol;

% Nonlinear parameter dependent problem (CHEBCOLLOC1)
pref.discretization = @chebcolloc1;

% Natural setup
x = chebfun('x', [-1 1]);
N = chebop(@(x, u, a) (1-x.^2).*u + .1*diff(u,2) + a.*exp(u));
% N = chebop(@(x, u, a) x.*u + .001*diff(u,2) + a.*diff(u));
N.lbc = @(u, a) [u + a + 1 ; diff(u)];
N.rbc = @(u, a) u - 1;
N.init = [x ; -1];
u = mldivide(N, chebfun(0), pref);
res = N(x, u);
Nlbc = N.lbc(u{1},u{2});
Nrbc = N.rbc(u{1},u{2});
err = norm(res) + Nlbc{1}(-1) + Nlbc{2}(-1) + Nrbc(1);

pass = err < tol;

end
