% Test file for bndfun constructor.

function pass = test_constructor(pref)

% Get preferences:
if ( nargin < 1 )
    pref = bndfun.pref;
end
% Set the tolerance:
pref = chebtech.pref(pref);
tol = 100*pref.chebtech.eps;

pass = zeros(1, 8); % Pre-allocate pass matrix.

%%
% Test on a scalar-valued function for extrapolation:
pref.chebtech.extrapolate = 1;
f = @(x) sin(x);
dom = [-2 7];
g = bndfun(f, dom, [], [], pref);
x = get(g,'points');
pass(1) = norm(f(x) - g.onefun.values, inf) < tol;

% Test on a scalar-valued function for interpolation:
f = @(x) sin(x)./x;
g = bndfun(f, dom, [], [], pref);
pass(2) = abs(1 - feval(g,0)) < tol;

%%
% Test on an array-valued function for extrapolation:
pref.chebtech.extrapolate = 1;
f = @(x) [sin(x) cos(x) exp(x)];
g = bndfun(f, dom, [], [], pref);
x = get(g,'points');
pass(3) = norm(f(x) - g.onefun.values, inf) < 30*tol;

% Test on an array-valued function for interpolation:
f = @(x) [sin(x)./x sin(x-3)./(x-3)];
g = bndfun(f, dom, [], [], pref);
gv = [feval(g,0) feval(g,3)];
pass(4) = norm(ones(1,2) - [gv(1) gv(4)], inf) < tol;

%%
% Test on chebtech preferences:
pref.chebtech.n = 50;
f = @(x) [sin(x) cos(x) exp(x)];
g = bndfun(f, dom, [], [], pref);
x = get(g,'points');
pass(5) = ( norm(f(x) - g.onefun.values, inf) < tol && size(g,1) == 50 );

%%
% Some other tests:
% reset preferences to factory values:
pref = bndfun.pref;
pref = chebtech.pref(pref);

% This should fail with an error:
try
    f = @(x) x + NaN;
    bndfun(f, dom, [], [], pref);
    pass(6) = false;
catch ME
    pass(6) = strcmp(ME.message, 'Too many NaNs to handle.');
end

% As should this:
try
    f = @(x) x + Inf;
    bndfun(f, dom, [], [], pref);
    pass(7) = false;
catch ME
    pass(7) = strcmp(ME.message, 'Too many NaNs to handle.');
end

% Test that the extrapolation option avoids endpoint evaluations.
pref.chebtech.extrapolate = 1;
try 
    bndfun(@(x) [F(x) F(x)], dom, [], [], pref);
    pass(8) = true;
catch ME %#ok<NASGU>
    pass(8) = false;
end

    function y = F(x)
        if ( any(abs(x) == 1) )
            error('Extrapolate should prevent endpoint evaluation.');
        end
        y = sin(x);
    end

end
