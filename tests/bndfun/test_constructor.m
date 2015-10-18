% Test file for bndfun constructor.

function pass = test_constructor(pref)

% Get preferences:
if ( nargin < 1 )
    pref = chebfunpref();
end

% Set the domain
dom = [-2 7];

% Generate a few random points to use as test values.
seedRNG(6178);
x = diff(dom) * rand(100, 1) + dom(1);

%%

% Test on a scalar-valued function for interpolation:
f = @(x) sin(x)./x;
g = bndfun(f, struct('domain', dom), pref);
pass(1) = abs(1 - feval(g, 0)) < 10*max(get(g, 'vscale')*eps);

%%

% Test on an array-valued function for interpolation:
f = @(x) [sin(x)./x sin(x - 3)./(x - 3)];
g = bndfun(f, struct('domain', dom), pref);
gv = [feval(g, 0) feval(g, 3)];
pass(2) = norm(ones(1, 2) - [gv(1) gv(4)], inf) < ...
    10*max(get(g, 'vscale')*eps);

%%
% Some other tests:

% This should fail with an error:
try
    f = @(x) x + NaN;
    bndfun(f, struct('domain', dom), pref);
    pass(3) = false;
catch ME
    pass(3) = strcmp(ME.message, 'Too many NaNs/Infs to handle.');
end

% As should this:
try
    f = @(x) x + Inf;
    bndfun(f, struct('domain', dom), pref);
    pass(4) = false;
catch ME
    pass(4) = strcmp(ME.message, 'Too many NaNs/Infs to handle.');
end

% Test that the extrapolation option avoids endpoint evaluations.
pref.techPrefs.extrapolate = 1;
try 
    bndfun(@(x) [F(x) F(x)], struct('domain', dom), pref);
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

%% Test on singular function:

powl = -0.5;
powr = -1.6;
op = @(x) (x - dom(1)).^powl.*(x - dom(2)).^powr.*sin(x);
pref.blowup = true;
data.domain = dom;
data.exponents = [powl powr];
f = bndfun(op, data, pref);
vals_f = feval(f, x);
vals_exact = feval(op, x);
err = vals_f-vals_exact;
pass(6) = ( norm(err, inf) < 1e3*eps*norm(vals_exact, inf) );

end
