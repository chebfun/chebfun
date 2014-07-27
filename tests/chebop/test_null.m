function pass = test_null(pref)

if ( nargin == 0 )
    pref = cheboppref();
end

tol = pref.errTol;

% Test 1:
N = chebop(@(u) diff(u), [0, pi]);
V = null(N);
err(1) = norm(N(V), inf);           % Check residual.
pass(1) = err(1) < tol;
pass(2) = size(V, 2) == 1;
err(3,1) = abs(V'*V - 1);           % Check orthanormality.
pass(3) = err(3) < tol;

% Test 2:
N = chebop(@(x, u) 0.2*diff(u, 3) - sin(3*x).*diff(u));
N.rbc = 1;
V = null(N);
plot(N(V))
err(4) = norm(N(V), inf);           % Check residual.
pass(4) = err(4) < 1e5*tol;         
pass(5) = size(V, 2) == 2;
err(6) = norm(V'*V - eye(2), inf);  % Check orthanomrality.
pass(6) = err(6) < tol;
err(7) = norm(feval(V, 1), inf);    % Check BC is satisfied.
pass(7) = err(7) < tol;

end


