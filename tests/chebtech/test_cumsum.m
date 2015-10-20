% Test file for chebtech/cumsum.m

function pass = test_cumsum(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebtech.techPref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

for n = 1:2
    if ( n == 1 )
        testclass = chebtech1();
    else 
        testclass = chebtech2();
    end

  %%
  % Spot-check antiderivatives for a couple of functions.  We verify that the
  % chebtech antiderivatives match the true ones up to a constant by checking 
  % that the standard deviation of the difference between the two on a large 
  % random grid is small. We also check that feval(cumsum(f), -1) == 0 each 
  % time.
  
  f = testclass.make(@(x) exp(x) - 1, [],  pref);
  F = cumsum(f);
  F_ex = @(x) exp(x) - x;
  err = std(feval(F, x) - F_ex(x));
  tol = 20*vscale(F)*eps;
  pass(n, 1) = (err < tol) && (abs(feval(F, -1)) < tol);
  
  f = testclass.make(@(x) 1./(1 + x.^2), [], pref);
  F = cumsum(f);
  F_ex = @(x) atan(x);
  err = feval(F, x) - F_ex(x);
  tol = 10*vscale(F)*eps;
  pass(n, 2) = (std(err) < tol) && (abs(feval(F, -1)) < tol);
  
  f = testclass.make(@(x) cos(1e4*x), [], pref);
  F = cumsum(f);
  F_ex = @(x) sin(1e4*x)/1e4;
  err = feval(F, x) - F_ex(x);
  tol = 5e4*vscale(F)*eps;
  pass(n, 3) = (std(err) < tol) && (abs(feval(F, -1)) < tol);
  
  z = exp(2*pi*1i/6);
  f = testclass.make(@(t) sinh(t*z), [], pref);
  F = cumsum(f);
  F_ex = @(t) cosh(t*z)/z;
  err = feval(F, x) - F_ex(x);
  tol = 10*vscale(F)*eps;
  pass(n, 4) = (std(err) < tol) && (abs(feval(F, -1)) < tol);
  
  %%
  % Check that applying cumsum() and direct construction of the antiderivative
  % give the same results (up to a constant).
  
  f = testclass.make(@(x) sin(4*x).^2, [], pref);
  F = testclass.make(@(x) 0.5*x - 0.0625*sin(8*x), [], pref);
  G = cumsum(f);
  err = G - F;
  tol = 10*vscale(G)*eps;
  values = err.coeffs2vals(err.coeffs); 
  pass(n, 5) = (std(values) < tol) && (abs(feval(G, -1)) < tol);
  
  %%
  % Check that diff(cumsum(f)) == f and that cumsum(diff(f)) == f up to a 
  % constant.
  
  f = testclass.make(@(x) x.*(x - 1).*sin(x) + 1, [], pref);
  g = diff(cumsum(f));
  err = feval(f, x) - feval(g, x);
  tol = 10*vscale(g)*eps;
  pass(n, 6) = (norm(err, inf) < 100*tol);
  h = cumsum(diff(f));
  err = feval(f, x) - feval(h, x);
  tol = 10*vscale(h)*eps;
  pass(n, 7) = (std(err) < tol)  && (abs(feval(h, -1)) < tol);
  
  %%
  % Check operation for array-valued chebtech objects.
  
  f = testclass.make(@(x) [sin(x) x.^2 exp(1i*x)], [], pref);
  F_exact = testclass.make(@(x) [(-cos(x)) (x.^3/3) (exp(1i*x)/1i)], [], pref);
  F = cumsum(f);
  err = std(feval(F, x) - feval(F_exact, x));
  tol = 10*max(vscale(F)*eps);
  pass(n, 8) = (norm(err, inf) < tol)  && all(abs(feval(F, -1)) < tol);
  
end

end
