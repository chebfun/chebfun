% Test file for bndfun/mldivide.m

function pass = test_mldivide(pref)

% Get preferences.
if (nargin < 1)
    pref = chebfunpref();
end

% Set a tolerance.  (pref.eps does not matter here.)
tol = 10*eps;

% Set the domain
dom = [-2 7];

% Generate a few random points to use as test values.
seedRNG(6178);
x = diff(dom) * rand(100, 1) + dom(1);

%%
% Basic correctness checks.

% We get a known exact solution in this case.
f = bndfun(@(x) sin(x), struct('domain', dom), pref);
g = f \ f;
err = f - g*f;
pass(1) = abs(g - 1) < tol;
pass(2) = max(abs(feval(err, x))) < tol;

% Same here.
f = bndfun(@(x) [sin(x) cos(x)], struct('domain', dom), pref);
g = bndfun(@(x) sin(x + pi/4), struct('domain', dom), pref);
h = f \ g;
err = g - f*h;
pass(3) = max(abs(h - [1/sqrt(2) ; 1/sqrt(2)])) < tol;
pass(4) = max(abs(feval(err, x))) < tol;

% A known least-squares solution.
f = bndfun(@(x) [ones(size(x)) x x.^2 x.^3], struct('domain', dom), pref);
g = bndfun(@(x) x.^4 + x.^3 + x + 1, struct('domain', dom), pref);
x = f \ g;
pass(5) = max(abs(x - [2469/70 ; -163/7 ; -141/7 ; 11])) < 1000*tol;

%%
% Check error conditions.

% mldivide doesn't work between a BNDFUN and a non-BNDFUN.
try
    f = bndfun(@(x) [sin(x) cos(x) exp(x)], struct('domain', dom), pref);
    f \ 2; %#ok<VUNUS>
    pass(6) = 0;
catch ME
    pass(6) = strcmp(ME.identifier, ...
       'CHEBFUN:BNDFUN:mldivide:bndfunMldivideUnknown');
end

end
