function pass = test_manualNewton(pref)
% Test that we can solve a problem with a manual Newton iteration.

if ( nargin == 0 )
    pref = cheboppref();
end

% This test is taken from the Chebfun guide.

% Specify a CHEBOP:
pref.discretization = @chebcolloc2;
L = chebop(@(x,u) 0.001*diff(u,2) - u.^3);
L.lbc = 1;
L.rbc = -1;
u = chebfun('-x');
nrmdu = Inf;
while nrmdu > 1e-10
   r = L*u;
   J = linearize(L,u);
   du = -mldivide(J,r,pref);
   u = u + du; 
   nrmdu = norm(du);
end

% If we get here, we must have converged:
pass = 1;

end
