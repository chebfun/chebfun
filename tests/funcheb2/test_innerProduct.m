% Test file for funcheb2/innerProduct.m

function pass = test_innerProduct(pref)

% Get preferences.
if ( nargin < 1 )
    pref = funcheb2.pref();
end

% Set a tolerance.  (pref.eps doesn't matter here.)
tol = 10*eps;

% Random numbers to use as arbitrary multiplicative constants.
alpha = -0.194758928283640 + 0.075474485412665i;
beta = -0.526634844879922 - 0.685484380523668i;

%%
% Spot-check a few known results.

f = funcheb2(@(x) sin(2*pi*x), pref);
g = funcheb2(@(x) cos(2*pi*x), pref);
pass(1) = abs(innerProduct(f, g)) < 10*max(f.epslevel, g.epslevel);

g = funcheb2(@(x) cos(4*pi*x), pref);
pass(2) = abs(innerProduct(f, g)) < 10*max(f.epslevel, g.epslevel);

f = funcheb2(@(x) exp(x), pref);
g = funcheb2(@(x) exp(-x), pref);
pass(3) = abs(innerProduct(f, g) - 2) < 10*max(f.epslevel, g.epslevel);

g = funcheb2(@(x) sin(x), pref);
exact = exp(1)*(sin(1) - cos(1))/2 - exp(-1)*(sin(-1) - cos(-1))/2;
pass(4) = abs(innerProduct(f, g) - exact) < 10*max(f.epslevel, g.epslevel);

%%
% Check a few known properties.

f = funcheb2(@(x) exp(x) - 1);
g = funcheb2(@(x) 1./(1 + 1i*x.^2));
h = funcheb2(@(x) sinh(x*exp(pi*1i/6)));

ip1 = innerProduct(alpha*f, beta*g);
ip2 = conj(alpha)*beta*innerProduct(f, g);
pass(5) = abs(ip1 - ip2) < tol;

ip1 = innerProduct(g, h);
ip2 = innerProduct(h, g);
pass(6) = abs(ip1 - conj(ip2)) < tol;

ip1 = innerProduct(f + g, h);
ip2 = innerProduct(f, h) + innerProduct(g, h);
pass(7) = abs(ip1 - ip2) < tol;

ip1 = innerProduct(f, g + h);
ip2 = innerProduct(f, g) + innerProduct(f, h);
pass(8) = abs(ip1 - ip2) < tol;

nf2 = innerProduct(f, f);
ng2 = innerProduct(g, g);
nh2 = innerProduct(h, h);
n2vals = [nf2 ; ng2 ; nh2];
pass(9) = isreal(n2vals) && all(n2vals >= 0);

%% 
% Check operation for vectorized funcheb2 objects.

f = funcheb2(@(x) [sin(x) cos(x)]);
g = funcheb2(@(x) [exp(x) 1./(1 + x.^2) airy(x)]);
ip = innerProduct(f, g);
exact = [0.663493666631241 0                 -0.135033172317858;
         1.933421496200713 1.365866063614065  0.592109441404267];
pass(10) = norm(ip(:) - exact(:), 'inf') < 10*max(f.epslevel, g.epslevel);

%%
% Check error conditions.
try
    ip = innerProduct(f, 2);
    pass(11) = false;
catch ME
    pass(11) = strcmp(ME.identifier, 'CHEBFUN:FUNCHEB2:InnerProduct:input');
end

% Can't take the inner product of a funcheb2 and a non-funcheb2.

end
