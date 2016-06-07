% Test file for trigtech/roots.m

function pass = test_roots(pref)

if ( nargin < 1 )
    pref = trigtech.techPref();
end

testclass = trigtech();

%% Simple test:
f = testclass.make(@(x) cos(5*pi*x), [], pref);
r = roots(f);
exact = (-0.9:0.2:0.9).';
pass(1) = norm(r-exact,Inf) < 1e1*length(f)*eps;

%% More complicated:
k = 20;
f = testclass.make(@(x) sin(sin(pi*k*x)), [], pref);
r = roots(f);
pass(2) = norm(r-(-k:k)'/k, inf) < length(f)*eps;

%% No real roots.
f = testclass.make( @(x) 3./(5 - 4*cos(3*pi*x)), ...
    [], pref);
r = roots(f);
pass(3) = isempty(r);

%% Test that complex roots are now allowed.
f = testclass.make(@(x) 2 + cos(pi*x), [], pref);
r = roots(f, 'complex', 1);
pass(4) = norm( feval(f,r) ) < vscale(f)*eps;

f = testclass.make(@(x) sin(100*pi*x));
r1 = roots(f, 'complex', 1);
r2 = roots(f);
pass(5) = numel(r1) == 200 & numel(r2) >= 201;

%% Test an array-valued function:
f = testclass.make(@(x) [sin(pi*x), cos(pi*x)], [], pref);
r = roots(f);
r2 = [-1 0 1 -.5 .5 NaN].';
pass(6) = all( r(:) - r2 < 10*length(f)*eps | isnan(r2) );

f = testclass.make(@(x) [cos(2*pi*x), sin(pi*x)], [], pref);
r = roots(f, 'complex', 1);
r = r(:);
[temp, id] = sort(real(r));
r = r(id);
r2 = [-0.75 -0.25 0 0.25 0.75 1 NaN NaN].';
pass(7) = all( abs(r(:) - sort(r2)) < 1e1*length(f)*eps | isnan(r2) );

end
