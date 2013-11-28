% Test file for chebfun constructor for integration with singfun.

function pass = test_constructor_singfun(pref)

if ( nargin == 0 )
    pref = chebpref();
end

%% specify the singularity by passing the exact exponents to preferences:
% No singularity detection involved: 

% define the domain:
dom = [-2 7];

% Generate a few random points to use as test values:
seedRNG(6178);
x = diff(dom) * rand(100, 1) + dom(1);

% create function handle:
pow = -0.5;
op = @(x) (x - dom(1)).^pow.*sin(x);

% specify the singularity in preference:
pref.singPrefs.exponents = [pow 0];

% construction:
f = chebfun(op, dom, pref);

% check values:
fval = feval(f, x);
vals_exact = feval(op, x);
err = fval - vals_exact;
pass(1) = ( norm(err, inf) < 1e3*get(f,'epslevel')*norm(vals_exact, inf) );

%% specify the singularity using Chebfun v4 syntax:
% No singularity detection involved: 

% define the domain:
dom = [-2 7];

% Generate a few random points to use as test values:
seedRNG(6178);
x = diff(dom) * rand(100, 1) + dom(1);

% create function handle:
pow = -0.5;
op = @(x) (x - dom(1)).^pow.*sin(x);

% construction:
f = chebfun(op, dom, 'exps', [pow 0]);

% check values:
fval = feval(f, x);
vals_exact = feval(op, x);
err = fval - vals_exact;
pass(2) = ( norm(err, inf) < 1e3*get(f,'epslevel')*norm(vals_exact, inf) );

%% specify the singularity by specifying the type of singularities to preferences:
% Singularity detection is involved: 

% define the domain:
dom = [-1 1];

% Generate a few random points to use as test values:
seedRNG(6178);
x = diff(dom) * rand(100, 1) + dom(1);

% create function handle:
pow = -0.5;
op = @(x) (x - dom(1)).^pow.*sin(x);

% specify the singularity in preference:
pref = chebpref();
pref.singPrefs.singType = {'sing', 'none'};

% construction:
f = chebfun(op, dom, pref);

% check values:
fval = feval(f, x);
vals_exact = feval(op, x);
err = fval - vals_exact;
pass(3) = ( norm(err, inf) < 1e3*get(f,'epslevel')*norm(vals_exact, inf) );

%% specify the singularity types using Chebfun v4 syntax:
% Singularity detection is involved: 

% define the domain:
dom = [-1 1];

% Generate a few random points to use as test values:
seedRNG(6178);
x = diff(dom) * rand(100, 1) + dom(1);

% create function handle:
pow = -0.5;
op = @(x) (x - dom(1)).^pow.*sin(x);

% construction:
f = chebfun(op, dom, 'blowup', 2);

% check values:
fval = feval(f, x);
vals_exact = feval(op, x);
err = fval - vals_exact;
pass(4) = ( norm(err, inf) < 1e3*get(f,'epslevel')*norm(vals_exact, inf) );

%% specify the singularity by specifying the type of singularities to preferences:
% Singularity detection is involved: 

% define the domain:
dom = [-1 1];

% Generate a few random points to use as test values:
seedRNG(6178);
x = diff(dom) * rand(100, 1) + dom(1);

% create function handle:
pow = -1;
op = @(x) (x - dom(1)).^pow.*sin(x);

% specify the singularity in preference:
pref = chebpref();
pref.singPrefs.singType = {'pole', 'none'};

% construction:
f = chebfun(op, dom, pref);

% check values:
fval = feval(f, x);
vals_exact = feval(op, x);
err = fval - vals_exact;
pass(5) = ( norm(err, inf) < 1e3*get(f,'epslevel')*norm(vals_exact, inf) );

%% specify the singularity types using Chebfun v4 syntax:
% Singularity detection is involved: 

% define the domain:
dom = [-1 1];

% Generate a few random points to use as test values:
seedRNG(6178);
x = diff(dom) * rand(100, 1) + dom(1);

% create function handle:
pow = -2;
op = @(x) (x - dom(1)).^pow.*sin(x);

% construction:
f = chebfun(op, dom, 'blowup', 1);

% check values:
fval = feval(f, x);
vals_exact = feval(op, x);
err = fval - vals_exact;
pass(6) = ( norm(err, inf) < 1e3*get(f,'epslevel')*norm(vals_exact, inf) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
% Test for default domain without specifying domain in chebfun constructor 
% calling sequence:

%% specify the singularity by specifying the type of singularities to preferences:
% Singularity detection is involved: 

% define the domain:
dom = [-1 1];

% Generate a few random points to use as test values:
seedRNG(6178);
x = diff(dom) * rand(100, 1) + dom(1);

% create function handle:
pow = -0.5;
op = @(x) (x - dom(1)).^pow.*sin(x);

% specify the singularity in preference:
pref = chebpref();
pref.singPrefs.singType = {'sing', 'none'};

% construction:
f = chebfun(op, pref);

% check values:
fval = feval(f, x);
vals_exact = feval(op, x);
err = fval - vals_exact;
pass(7) = ( norm(err, inf) < 1e3*get(f,'epslevel')*norm(vals_exact, inf) );

%% specify the singularity types using Chebfun v4 syntax:
% Singularity detection is involved: 

% define the domain:
dom = [-1 1];

% Generate a few random points to use as test values:
seedRNG(6178);
x = diff(dom) * rand(100, 1) + dom(1);

% create function handle:
pow = -0.5;
op = @(x) (x - dom(1)).^pow.*sin(x);

% construction:
f = chebfun(op, 'blowup', 2);

% check values:
fval = feval(f, x);
vals_exact = feval(op, x);
err = fval - vals_exact;
pass(8) = ( norm(err, inf) < 1e3*get(f,'epslevel')*norm(vals_exact, inf) );

%% specify the singularity by specifying the type of singularities to preferences:
% Singularity detection is involved: 

% define the domain:
dom = [-1 1];

% Generate a few random points to use as test values:
seedRNG(6178);
x = diff(dom) * rand(100, 1) + dom(1);

% create function handle:
pow = -1;
op = @(x) (x - dom(1)).^pow.*sin(x);

% specify the singularity in preference:
pref = chebpref();
pref.singPrefs.singType = {'pole', 'none'};

% construction:
f = chebfun(op, pref);

% check values:
fval = feval(f, x);
vals_exact = feval(op, x);
err = fval - vals_exact;
pass(9) = ( norm(err, inf) < 1e3*get(f,'epslevel')*norm(vals_exact, inf) );

%% specify the singularity types using Chebfun v4 syntax:
% Singularity detection is involved: 

% define the domain:
dom = [-1 1];

% Generate a few random points to use as test values:
seedRNG(6178);
x = diff(dom) * rand(100, 1) + dom(1);

% create function handle:
pow = -2;
op = @(x) (x - dom(1)).^pow.*sin(x);

% construction:
f = chebfun(op, 'blowup', 1);

% check values:
fval = feval(f, x);
vals_exact = feval(op, x);
err = fval - vals_exact;
pass(10) = ( norm(err, inf) < 1e3*get(f,'epslevel')*norm(vals_exact, inf) );

end