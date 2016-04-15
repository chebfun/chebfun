% Test file for chebfun/innerProduct.m

function pass = test_innerProduct(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

% Set a domain
dom = [-2 7];

%% Test for singular function: piecewise smooth chebfun - splitting on.

pow1 = -0.3;
pow2 = -0.5;
op1 = @(x) (x - dom(2)).^pow1.*sin(100*x);
op2 = @(x) (x - dom(2)).^pow2.*cos(300*x);
f = chebfun(op1, dom, 'exps', [0 pow1], 'splitting', 'on');
g = chebfun(op2, dom, 'exps', [0 pow2], 'splitting', 'on');
I = innerProduct(f,g);
I_exact = 0.35838148154346034 - 0.26037938759089226i;
pass(1) = ( abs(I-I_exact) < 1e5*eps*...
    abs(I_exact) );


%% Test for functions defined on unbounded domains:

%% Functions on [-inf inf]:

% Set the domain:
dom = [-Inf Inf];

opf = @(x) 2-exp(-x.^2);
opg = @(x) exp(-x.^2);
f = chebfun(opf, dom);
g = chebfun(opg, dom);

I = innerProduct(f, g);
IExact = (sqrt(pi)*(4 - sqrt(2)))/2;
err = abs(I - IExact);
pass(2) = err < 2e7*max(eps*get(f,'vscale'), ...
    eps*get(g,'vscale'));

%% Functions on [a inf]:

% Set the domain:
dom = [1 Inf];

opf = @(x) x;
opg = @(x) exp(-x);
f = chebfun(opf, dom, 'exps', [0 1]);
g = chebfun(opg, dom);
warning('off', 'CHEBFUN:UNBNDFUN:sum:slowDecay');
I = innerProduct(f, g);
warning('on', 'CHEBFUN:UNBNDFUN:sum:slowDecay');
IExact = 2*exp(-1);
err = abs(I - IExact);
pass(3) = err < 1e9*max(eps*get(f,'vscale'), ...
    eps*get(g,'vscale'));

%% Functions on [-inf b]:

% Set the domain:
dom = [-Inf -3*pi];

opf = @(x) 1./x;
opg = @(x) 2./x;
f = chebfun(opf, dom);
g = chebfun(opg, dom);
I = innerProduct(f, g);
IExact = 2/(3*pi);
err = abs(I - IExact);
pass(4) = err < 1e5*eps*get(f,'vscale');

%% From #559
x = chebfun('x', [-1 1], pref);
X = [x x];
S = sqrt(abs(x));
I = innerProduct(X, S);
err = norm(I, inf);
tol = eps*get(x,'vscale');
pass(5) = err < tol;
E = exp(X);
I = innerProduct(E, S);
err = norm(I - 1.634574774192848, inf);
tol = max(eps*get(E,'vscale'));
pass(6) = err < tol;

% Test inner products between chebfuns based on mixed techs.  (See #1784.)
f = chebfun(@(x) sin(pi*x));
g = chebfun(@(x) sin(pi*x), 'trig');
err = abs(innerProduct(f, g) - 1);
tol = 10*eps*get(f, 'vscale');
pass(7) = err < tol;

f = chebfun(@(x) [sin(pi*x) cos(pi*x)]);
g = chebfun(@(x) [sin(pi*x) cos(pi*x)], 'trig');
err = norm(innerProduct(f, g) - [1 0 ; 0 1], Inf);
tol = 10*eps*get(f, 'vscale');
pass(8) = err < tol;

f = cheb2quasi(f);
g = cheb2quasi(g);
err = norm(innerProduct(f, g) - [1 0 ; 0 1], Inf);
tol = 10*eps*get(f, 'vscale');
pass(9) = err < tol;

end
