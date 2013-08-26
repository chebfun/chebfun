% Test file for chebtech/cumsum.m

% [TODO]:  Fix test tolerances once we've figured out how to properly update
% epslevels in chebtech.cumsum().

function pass = test_cumsum(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebtech.pref;
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

pass = zeros(2, 10); % Pre-allocate pass matrix
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
  
  f = testclass.make(@(x) exp(x) - 1, [], [],  pref);
  F = cumsum(f);
  F_ex = @(x) exp(x) - x;
  err = feval(F, x) - F_ex(x);
  tol = 10*F.vscale.*F.epslevel;
  pass(n, 1) = (std(err) < tol) && (abs(feval(F, -1)) < tol);
  
  f = testclass.make(@(x) 1./(1 + x.^2), [], [], pref);
  F = cumsum(f);
  F_ex = @(x) atan(x);
  err = feval(F, x) - F_ex(x);
  tol = 10*F.vscale.*F.epslevel;
  pass(n, 2) = (std(err) < tol) && (abs(feval(F, -1)) < tol);
  
  f = testclass.make(@(x) cos(1e4*x), [], [], pref);
  F = cumsum(f);
  F_ex = @(x) sin(1e4*x)/1e4;
  err = feval(F, x) - F_ex(x);
  tol = 1e4*F.vscale.*F.epslevel;
  pass(n, 3) = (std(err) < tol) && (abs(feval(F, -1)) < tol);
  
  z = exp(2*pi*1i/6);
  f = testclass.make(@(t) sinh(t*z), [], [], pref);
  F = cumsum(f);
  F_ex = @(t) cosh(t*z)/z;
  err = feval(F, x) - F_ex(x);
  tol = 10*F.vscale.*F.epslevel;
  pass(n, 4) = (std(err) < tol) && (abs(feval(F, -1)) < tol);
  
  %%
  % Check that applying cumsum() and direct construction of the antiderivative
  % give the same results (up to a constant).
  
  f = testclass.make(@(x) sin(4*x).^2, [], [], pref);
  F = testclass.make(@(x) 0.5*x - 0.0625*sin(8*x), [], [], pref);
  G = cumsum(f);
  err = G - F;
  tol = 10*G.vscale.*G.epslevel;
  pass(n, 5) = (std(err.values) < tol) && (abs(feval(G, -1)) < tol);
  
  %%
  % Check that diff(cumsum(f)) == f and that cumsum(diff(f)) == f up to a 
  % constant.
  
  f = testclass.make(@(x) x.*(x - 1).*sin(x) + 1, [], [], pref);
  g = diff(cumsum(f));
  err = feval(f, x) - feval(g, x);
  tol = 10*g.vscale.*g.epslevel;
  pass(n, 6) = (norm(err, inf) < 100*tol);
  h = cumsum(diff(f));
  err = feval(f, x) - feval(h, x);
  tol = 10*h.vscale.*h.epslevel;
  pass(n, 7) = (std(err) < tol)  && (abs(feval(h, -1)) < tol);
  
  %%
  % Check operation for array-valued chebtech objects.
  
  f = testclass.make(@(x) [sin(x) x.^2 exp(1i*x)], [], [], pref);
  F_exact = testclass.make(@(x) [(-cos(x)) (x.^3/3) (exp(1i*x)/1i)], [], [], pref);
  F = cumsum(f);
  err = std(feval(F, x) - feval(F_exact, x));
  tol = 10*max(F.vscale.*F.epslevel);
  pass(n, 8) = (norm(err, inf) < tol)  && all(abs(feval(F, -1)) < tol);
  
  %%
  % Check operation for second and third order cumsums.
  f = testclass.make(@(x) sin(x), [], [], pref);
  F2_exact = testclass.make(@(x) -sin(x)+x.*cos(-1)+sin(-1)+cos(-1), ...
      [], [], pref);
  F2 = cumsum(f,2);
  err = std(feval(F2, x) - feval(F2_exact, x));
  tol = 10*F2.vscale.*F2.epslevel;
  pass(n, 9) = (norm(err, inf) < tol)  && abs(feval(F2, -1) < tol);
  
  F3_exact = testclass.make(@(x) cos(x)+x.^2*cos(-1)/2+x*(sin(-1)+cos(-1)) - ...
      (cos(-1)/2-sin(-1)), [], [], pref);
  F3 = cumsum(f,3);
  err = std(feval(F3, x) - feval(F3_exact, x));
  tol = 10*F3.vscale.*F3.epslevel;
  pass(n, 10) = (norm(err, inf) < tol)  && abs(feval(F3, -1) < tol);
  
  
end

end
