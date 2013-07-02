% Test file for bndfun constructor.

function pass = test_constructor(pref)

% Get preferences:
if ( nargin < 1 )
    pref = fun.pref;
end

% Set the tolerance:
tol = 10*pref.fun.eps;

% Set the domain
dom = [-2 7];

pass = zeros(1, 5); % Pre-allocate pass matrix.

%%

% Test on a scalar-valued function for interpolation:
f = @(x) sin(x)./x;
g = bndfun(f, dom, [], [], pref);
pass(1) = abs(1 - feval(g, 0)) < max(g.onefun.vscale)*tol;

%%

% Test on an array-valued function for interpolation:
f = @(x) [sin(x)./x sin(x - 3)./(x - 3)];
g = bndfun(f, dom, [], [], pref);
gv = [feval(g, 0) feval(g, 3)];
pass(2) = norm(ones(1, 2) - [gv(1) gv(4)], inf) < max(g.onefun.vscale)*tol;


%%
% Some other tests:
% Reset preferences to factory values:
pref = bndfun.pref;
pref = chebtech.pref(pref);

% This should fail with an error:
try
    f = @(x) x + NaN;
    bndfun(f, dom, [], [], pref);
    pass(3) = false;
catch ME
    pass(3) = strcmp(ME.message, 'Too many NaNs/Infs to handle.');
end

% As should this:
try
    f = @(x) x + Inf;
    bndfun(f, dom, [], [], pref);
    pass(4) = false;
catch ME
    pass(4) = strcmp(ME.message, 'Too many NaNs/Infs to handle.');
end

% Test that the extrapolation option avoids endpoint evaluations.
pref.chebtech.extrapolate = 1;
try 
    bndfun(@(x) [F(x) F(x)], dom, [], [], pref);
    pass(5) = true;
catch ME %#ok<NASGU>
    pass(5) = false;
end

    function y = F(x)
        if ( any(abs(x) == 1) )
            error('Extrapolate should prevent endpoint evaluation.');
        end
        y = sin(x);
    end

end