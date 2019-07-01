function pass = test_eigs_piecewise(pref)

% This test is taken from #1074.

if ( nargin == 0 )
    pref = cheboppref();
end

x = chebfun('x');
ep = 0.25;
F = (abs(x) < ep)/(2*ep);
L = chebop(@(x,u) diff(u,2), [-1 1], 0, 0);
M = chebop(@(x,u) F(x).*u);

pref.discretization = @chebcolloc2;
ec2 = eigs(L, M, pref);

pref.discretization = @chebcolloc1;
ec1 = eigs(L, M, pref);

pref.discretization = @ultraS;
eus = eigs(L, M, pref);

refSoln =    1.0e+02 * ...
 [  -4.987961490158024
  -3.211336584876317
  -1.829388534636815
  -0.841876384196023
  -0.247291299790162
  -0.023950791540263];

tol = 1e4*pref.bvpTol;
err(1) = norm(ec2 - refSoln, inf);
err(2) = norm(ec1 - refSoln, inf);
err(3) = norm(eus - refSoln, inf);

pass = err < tol;

end

