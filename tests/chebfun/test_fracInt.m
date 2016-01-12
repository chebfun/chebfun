% Test file for fractional calculus.

function pass = test_fracInt(pref)

if ( nargin == 0 ) 
    pref = chebfunpref();
end

tol = 100*pref.eps;
dom = [0, 1];

% Test values:
xx = linspace(dom(1)+.1, dom(2)-.1, 10)';

%% Polynomials
x = chebfun('x', dom);
q = sqrt(2)/2;
k = 1;
for n = [1, 4]
    
    U = diff(x.^n, q);
    tru = gamma(n+1)./gamma(n+1-q)*chebfun(@(x) x.^(n-q), dom, 'exps', [n-q, 0]); 
    
    err(k) = norm(tru(xx) - U(xx), inf);
    tol(k) = 10*eps*vscale(U)*hscale(U);
    k = k + 1;
    
end

%% Exponential
u = chebfun('exp(x)', dom);
trueC = chebfun('erf(sqrt(x)).*exp(x) + 1./sqrt(pi*x)', dom, 'exps', [-.5 0]);
trueRL = chebfun('erf(sqrt(x)).*exp(x)', dom, 'exps', [.5, 0]);

% RL
U = diff(u, .5, 'Caputo');
err(3) = norm(trueRL(xx) - U(xx), inf);
tol(3) = 10*eps*vscale(U)*hscale(U);

% Caputo
U = diff(u, .5, 'RL');
err(4) = norm(trueC(xx) - U(xx), inf);
tol(4) = 1e2*eps*vscale(U)*hscale(U);
    

%% Integrate twice:
xx = linspace(-sqrt(2)*pi+.1, pi-.1, 10).';
f = chebfun(@sin, [-sqrt(2)*pi, pi]);
F = cumsum(f);
G = fracInt(fracInt(f, .3), .7);
err(5) = norm(feval(F, xx) - feval(G, xx), inf);
tol(5) = 10*eps*vscale(G)*hscale(G);

%% Differentiate twice:
xx = linspace(-sqrt(2)*pi+.1, pi-.1, 10)';
f = chebfun(@sin, [-sqrt(2)*pi, pi]);
F = diff(f);
G = fracDiff(fracDiff(f, .3), .7);
err(6) = norm(feval(F, xx) - feval(G, xx), inf);
tol(6) = 1e2*eps*vscale(G)*hscale(G);


%%
pass = err < tol;

end
