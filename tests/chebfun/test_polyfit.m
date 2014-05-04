% Test file for POLYFIT.

function pass = test_polyfit(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

seedRNG(7681);

%% Test CONTINUOUSPOLYFIT on [-1 1]:

% Check the 'nothing to do' case:
F = chebfun(@(x) x.^2 + 3*x.^4, [-1  1]);
p = polyfit(F, 5);
err = norm(p - F);
pass(1) = err < 10*epslevel(p)*vscale(p);

% Check the 'nothing to do' case with breakpoints:
F = chebfun(@(x) x.^2 + 3*x.^4, [-1 -0.5 0 0.5 1]);
p = polyfit(F, 6);
err = norm(p - F);
pass(2) = err < 10*epslevel(p)*vscale(p);

% Fit a quartic chebfun:
F = chebfun(@(x) x.^2 + 3*x.^4, [-1 1]);
p = polyfit(F, 4);
err = norm(p - F);
pass(3) = err < 10*epslevel(p)*vscale(p);

% Fit a quartic chebfun with breakpoints:
F = chebfun(@(x) x.^2 + 3*x.^4, [-1 -0.5 0 0.5 1]);
p = polyfit(F, 5);
err = norm(p - F);
pass(4) = err < 10*epslevel(p)*vscale(p);

% Fit a row quartic chebfun:
F = chebfun(@(x) x.^2 + 3*x.^4, [-1 1]).';
p = polyfit(F, 4);
err = norm(p - F);
pass(5) = err < 10*epslevel(p)*vscale(p);

% Fit a row quartic chebfun with breakpoints:
F = chebfun(@(x) x.^2 + 3*x.^4, [-1 -0.5 0 0.5 1]).';
p = polyfit(F, 4);
err = norm(p - F);
pass(6) = err < 10*epslevel(p)*vscale(p);

%% Test CONTINUOUSPOLYFIT on [0 500]:

% Fit a quartic chebfun with breakpoint:
F = chebfun(@(x) 3*x.^2, [0 100 500]);
p = polyfit(F, 2);
err = norm(p - F);
pass(7) = err < 1e2*epslevel(p)*vscale(p);

%% Test CONTINUOUSPOLYFIT for array-valued inputs:
F = chebfun(@(x) [3*x.^3 (x.^2 + 3*x.^4)], [-1 -0.5 0 0.5 1]);
p = polyfit(F, 4);
err = norm(p - F);
pass(8) = err < 10*epslevel(p)*vscale(p);

end
