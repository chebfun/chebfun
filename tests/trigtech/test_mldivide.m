% Test file for trigtech/mldivide.m

function pass = test_mldivide(pref)

% Get preferences.
if ( nargin < 1 )
    pref = trigtech.techPref();
end

testclass = trigtech();

%%
% Basic correctness checks.

% We get a known exact solution in this case.
f = testclass.make(@(x) cos(sin(pi*x)), [], pref);
x = f \ f;
err = f - x*f;
pass(1) = abs(x - 1) < 10*f.vscale.*f.epslevel;
pass(2) = max(abs(err.coeffs(:))) < 10*f.vscale.*f.epslevel;

% Same here.
f = testclass.make(@(x) [sin(pi*x) cos(pi*x)], [], pref);
g = testclass.make(@(x) sin(pi*x + pi/4), [], pref);
tol_f = 10*max(f.vscale.*f.epslevel);
tol_g = 10*max(f.vscale.*f.epslevel);
x = f \ g;
err = g - f*x;
pass(3) = max(abs(x - [1/sqrt(2) ; 1/sqrt(2)])) < max(tol_f, tol_g);
pass(4) = max(abs(err.values(:))) < max(tol_f, tol_g);

% A known least-squares solution.
f = testclass.make(@(x) [ones(size(x)) cos(pi*x) sin(pi*x)], [], pref);
g = testclass.make(@(x) cos(pi*x), [], pref);
tol_f = 10*max(f.vscale.*f.epslevel);
tol_g = 10*max(f.vscale.*f.epslevel);
x = f \ g;
pass(5) = max(abs(x - [0 ; 1 ; 0])) < max(tol_f, tol_g);

%%
% Check error conditions.

% MLDIVIDE doesn't work between a TRIGTECH and a non-TRIGTECH.
try
    f = testclass.make(@(x) [sin(pi*x) cos(pi*x) exp(1i*pi*x)], [], pref);
    f \ 2; %#ok<VUNUS>
    pass(6) = 0;
catch ME
    pass(6) = strcmp(ME.identifier, ...
        'CHEBFUN:TRIGTECH:mldivide:trigtechMldivideUnknown');
end

end
