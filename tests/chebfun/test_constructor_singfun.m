% Test file for chebfun constructor for singular function.

function pass = test_constructor_singfun(pref)

if ( nargin == 0 )
    pref = chebfunpref();
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

% construction:
f = chebfun(op, dom, 'exps', [pow 0]);

% check values:
fval = feval(f, x);
vals_exact = feval(op, x);
err = fval - vals_exact;
pass(1) = ( norm(err, inf) < 1e3*eps*norm(vals_exact, inf) );

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
pass(2) = ( norm(err, inf) < 1e3*eps*norm(vals_exact, inf) );

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

% construction:
f = chebfun(op, dom, 'blowup', 2, ...
    'singType', {'sing', 'none'});

% check values:
fval = feval(f, x);
vals_exact = feval(op, x);
err = fval - vals_exact;
pass(3) = ( norm(err, inf) < 1e3*eps*norm(vals_exact, inf) );

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
pass(4) = ( norm(err, inf) < 1e3*eps*norm(vals_exact, inf) );

% Same thing but using 'blowup' and flag '2' instead:
f = chebfun(op, dom, 'blowup', 2);
fval = feval(f, x);
vals_exact = feval(op, x);
err = fval - vals_exact;
pass(5) = ( norm(err, inf) < 1e3*eps*norm(vals_exact, inf) );

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

% construction:
f = chebfun(op, dom, 'blowup', 2, ...
    'singType', {'pole', 'none'});

% check values:
fval = feval(f, x);
vals_exact = feval(op, x);
err = fval - vals_exact;
pass(6) = ( norm(err, inf) < 1e3*eps*norm(vals_exact, inf) );

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
pass(7) = ( norm(err, inf) < 1e3*eps*norm(vals_exact, inf) );

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

% construction:
f = chebfun(op, 'blowup', 2, ...
    'singType', {'sing', 'none'});

% check values:
fval = feval(f, x);
vals_exact = feval(op, x);
err = fval - vals_exact;
pass(8) = ( norm(err, inf) < 1e3*eps*norm(vals_exact, inf) );

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
pass(9) = ( norm(err, inf) < 1e3*eps*norm(vals_exact, inf) );

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

% construction:
f = chebfun(op, 'blowup', 2, ...
    'singType', {'pole', 'none'});

% check values:
fval = feval(f, x);
vals_exact = feval(op, x);
err = fval - vals_exact;
pass(10) = ( norm(err, inf) < 1e3*eps*norm(vals_exact, inf) );

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
pass(11) = ( norm(err, inf) < 1e3*eps*norm(vals_exact, inf) );

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
    result(j) = ( norm(err, inf) < 1e3*eps*norm(vals_exact, inf) );
end

pass(12) = all( result );

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
pass(13) = all( norm(err, inf) < 1e4*eps*norm(vals_exact, inf) );


%% Multipole subdomains with blowup flag set 1:

% Set the domain:
dom = [-2 0 7];
op1 = @(x) sin(x)./((x-dom(1)).*(dom(2)-x));
op2 = @(x) cos(x)./((x-dom(2)).^2.*(dom(3)-x));
op = {op1 op2};
f = chebfun(op, dom, 'blowup', 1);

% check values:

% Generate a few random points to use as test values:
x1 = diff(dom(1:2)) * rand(100, 1) + dom(1);
x2 = diff(dom(2:3)) * rand(100, 1) + dom(2);

fval1 = feval(f, x1);
vals_exact1 = feval(op1, x1);
err1 = fval1 - vals_exact1;
fval2 = feval(f, x2);
vals_exact2 = feval(op2, x2);
err2 = fval2 - vals_exact2;

pass(14) = norm([err1; err2], inf) < 1e2*eps* ...
    norm([vals_exact1; vals_exact2], inf);


% Exponents is set to NaN:
op = @(x) exp(x).*sqrt(1+x)./(1-x).^2;
f = chebfun(op, 'exps', [.5 NaN]);
x = 2 * rand(100, 1) - 1;
fx = feval(f, x);
f_exact = op(x);
err = fx - f_exact;
pass(15) = norm(err, inf) < 1e1*eps*norm(f_exact, inf);

%% See #1486

f = chebfun(@(x) tan(x),[0 pi],'splitting','on','blowup','on');
pass(16) = norm(f, inf) > 0;
f = chebfun(@(x) tan(x),[5*pi,6*pi], 'splitting', 'on', 'blowup', 1);
pass(17) = numel(f.domain) == 3;

end
