% Test file for trigtech/cumsum.m
function pass = test_cumsum(pref)

% Get preferences.
if ( nargin < 1 )
    pref = trigtech.techPref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

testclass = trigtech();

%%
% Spot-check antiderivatives for a couple of functions.  We verify that the
% trigtech antiderivatives match the true ones up to a constant by checking 
% that the standard deviation of the difference between the two on a large 
% random grid is small. We also check that feval(cumsum(f), -1) == 0 each 
% time.

% Note that trigtech's cumsum only works when mean of the function is zero.
% Thus in all these tests, the functions have this property.  There is one
% test at the end that verifies an error message is given when cumsum is
% applied to a function without zero mean.
 
k = 2; a = 0;
f = testclass.make(@(x) sin(k*pi*(x-a)).*cos((k+1)*pi*(x-a)), [],  pref);
F = cumsum(f);
F_ex = @(x) ((2*k+1)*cos(pi*(x-a))-cos((2*k+1)*pi*(x-a)))/(2*(pi+2*k*pi));
err = feval(F, x) - F_ex(x);
tol = 10*F.vscale.*F.epslevel;
pass(1) = (std(err) < tol) && (abs(feval(F, -1)) < tol);

k = 200; a = 0.17;
f = testclass.make(@(x) sin(k*pi*(x-a)).*cos((k+1)*pi*(x-a)), [],  pref);
F = cumsum(f);
F_ex = @(x) ((2*k+1)*cos(pi*(x-a))-cos((2*k+1)*pi*(x-a)))/(2*(pi+2*k*pi));
err = feval(F, x) - F_ex(x);
tol = 100*F.vscale.*F.epslevel;
pass(2) = (std(err) < tol) && (abs(feval(F, -1)) < tol);

k1 = 5; a1 = -0.33; k2 = 40; a2 = 0.17;
f = testclass.make(@(x) sin(k1*pi*(x-a1)).*cos((k1+1)*pi*(x-a1)) + 1i*sin(k2*pi*(x-a2)).*cos((k2+1)*pi*(x-a2)), [],  pref);
F = cumsum(f);
F_ex = @(x) ((2*k1+1)*cos(pi*(x-a1))-cos((2*k1+1)*pi*(x-a1)))/(2*(pi+2*k1*pi)) + 1i*((2*k2+1)*cos(pi*(x-a2))-cos((2*k2+1)*pi*(x-a2)))/(2*(pi+2*k2*pi));
err = feval(F, x) - F_ex(x);
tol = 100*F.vscale.*F.epslevel;
pass(3) = (std(err) < tol) && (abs(feval(F, -1)) < tol);

%%
% Check that diff(cumsum(f)) == f and that cumsum(diff(f)) == f up to a 
% constant.

f = testclass.make(@(x) sin(4*pi*cos(pi*x)), [], pref);
g = diff(cumsum(f));
err = feval(f, x) - feval(g, x);
tol = 10*g.vscale.*g.epslevel;
pass(4) = (norm(err, inf) < 100*tol);
h = cumsum(diff(f));
err = feval(f, x) - feval(h, x);
tol = 10*h.vscale.*h.epslevel;
pass(5) = (std(err) < tol)  && (abs(feval(h, -1)) < tol);

%%
% Check operation for array-valued trigtech objects.

f = testclass.make(@(x) [sin(4*pi*cos(2*pi*x)) sin(3*pi*x)], [], pref);
g = diff(cumsum(f));
err = feval(f, x) - feval(g, x);
tol = 10*g.vscale.*g.epslevel;
pass(6) = all(max(abs(err)) < 100*tol);
h = cumsum(diff(f));
err = feval(f, x) - feval(h, x);
tol = 10*h.vscale.*h.epslevel;
pass(7) = all((std(err) < tol))  && all(abs(feval(h, -1)) < tol);
  
%%
% Check that an error is thrown when the mean of the trigtech is not zero.

f = testclass.make(@(x) exp(cos(pi*x)), [], pref);
try
    g = cumsum(f);
    pass(8) = false;
catch ME
    pass(8) = strcmp(ME.identifier, 'CHEBFUN:TRIGTECH:cumsum:meanNotZero');
end

%%
% Check that an error is thrown when the mean just one of the means of 
% an array-valued trigtech is not zero.
f = testclass.make(@(x) [sin(4*pi*cos(pi*x)) exp(cos(pi*x))], [], pref);
try
    g = cumsum(f);
    pass(9) = false;
catch ME
    pass(9) = strcmp(ME.identifier, 'CHEBFUN:TRIGTECH:cumsum:meanNotZero');
end

end
