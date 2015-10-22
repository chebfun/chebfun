% Test file for trigtech/flipud.m

function pass = test_flipud(pref)

% Get preferences.
if ( nargin < 1 )
    pref = trigtech.techPref();
end
    
testclass = trigtech();

%%
% Conduct a few very straightforward tests.
f = testclass.make(@(x) sin(pi*x), [], pref);
g = testclass.make(@(x) -sin(pi*x), [], pref);
h = flipud(f);
pass(1) = norm(g.coeffs - h.coeffs, inf) < 10*vscale(h)*eps;

f = testclass.make(@(x) [sin(sin(pi*x)), exp(1i*pi*x)], [], pref);
g = testclass.make(@(x) [-sin(sin(pi*x)), exp(-1i*pi*x)], [], pref);
h = flipud(f);
pass(2) = norm(g.coeffs - h.coeffs, inf) < 100*max(vscale(h)*eps);

end
