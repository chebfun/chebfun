% Test file for chebfun constructor for singular function:.

function pass = test_constructor_singfun(pref)

if ( nargin == 0 )
    pref = chebpref();
end

%% Specify the singularity by passing the exact exponents to preferences:
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

%% Specify the singularity using Chebfun v4 syntax:
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

%% Specify the singularity by specifying the type of singularities to preferences:
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

%% Specify the singularity types using Chebfun v4 syntax:
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

%% Specify the singularity by naming the type of singularities to preferences:
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

%% Specify the singularity types using Chebfun v4 syntax:
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

%% Specify the singularity types using Chebfun v4 syntax:
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

%% Specify the singularity by specifying the type of singularities to preferences:
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

%% Specify the singularity types using Chebfun v4 syntax:
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

%% Piecewise smooth chebfun: a mix of functions with finite and infinite values.
% define the domain:
dom = [-2 -1 0 1];

op1 = @(x) sin(x);
op2 = @(x) 1./(1+x);
op3 = @(x) x+1;
op = {op1, op2, op3};
f = chebfun(op, dom, 'exps', [0 0 -1 0 0 0]);

% check values:
result = zeros(1,3);
for j = 1:3
    % define the domain:
    curDom = dom(j:j+1);
    
    % Generate a few random points to use as test values:
    x = diff(curDom) * rand(100, 1) + curDom(1);
    
    fval = feval(f, x);
    vals_exact = feval(op{j}, x);
    err = fval - vals_exact;
    result(j) = ( norm(err, inf) < 1e3*get(f,'epslevel')*norm(vals_exact, inf) );
end

pass(11) = all( result );

%% Piecewise smooth chebfun: splitting on.

% define the domain:
dom = [-2 3];

op = @(x) sin(300*x)./((x-dom(1)).*(x-dom(2)));
f = chebfun(op, dom, 'exps', [-1 -1], 'splitting', 'on');

% check values:

% Generate a few random points to use as test values:
x = diff(dom) * rand(100, 1) + dom(1);

fval = feval(f, x);
vals_exact = feval(op, x);
err = fval - vals_exact;
pass(12) = all( norm(err, inf) < 1e1*get(f,'epslevel')*norm(vals_exact, inf) );

end