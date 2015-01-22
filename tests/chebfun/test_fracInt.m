% Test file for fractional calculus.

function pass = test_fracInt(pref)

if ( nargin == 0 ) 
    pref = chebfunpref();
end

tol = 100*pref.eps;
dom = [0, 1];

% Test values:
xx = linspace(dom(1)+.1, dom(2)-.1, 100);

%% Polynomials
x = chebfun('x', dom);
q = sqrt(2)/2;
k = 1;
for n = [1, 4]
    
    xnpq = diff(x.^n, q, 'Caputo');
    tru = gamma(n+1)./gamma(n+1-q)*chebfun(@(x) x.^(n-q), dom, 'exps', [n-q, 0]); 
    
    err(k) = norm(tru(xx) - xnpq(xx), inf);
    k = k + 1;
    
end

%% Exponential
u = chebfun('exp(x)', dom);
trueC = chebfun('erf(sqrt(x)).*exp(x) + 1./sqrt(pi*x)', dom, 'exps', [-.5 0]);
trueRL = chebfun('erf(sqrt(x)).*exp(x)', dom, 'exps', [.5, 0]);

% RL
up05a = diff(u, .5, 'RL');
err(3) = norm(trueRL(xx) - up05a(xx), inf);

% Caputo
up05b = diff(u, .5, 'Caputo');
err(4) = norm(trueC(xx) - up05b(xx), inf);

%%

pass = err < tol;

end