% Test file for bndfun/innerProduct.m

function pass = test_innerProduct(pref)

% Get preferences.
if ( nargin < 1 )
    pref = bndfun.pref;
end

% Set a tolerance.  (pref.eps doesn't matter here.)
tol = 10*eps;

% Set a domain
dom = [-2 7];

% Fixed arbitrary numbers to use as multiplicative constants.
alpha = -0.194758928283640 + 0.075474485412665i;
beta = -0.526634844879922 - 0.685484380523668i;

pass = zeros(1, 11); % Pre-allocate pass matrix

%%
% Spot-check a few known results.

f = bndfun(@(x) sin(2*pi*x), dom, [], [], pref);
g = bndfun(@(x) cos(2*pi*x), dom, [], [], pref);
pass(1) = abs(innerProduct(f, g)) < 10*max(f.onefun.epslevel, g.onefun.epslevel);

g = bndfun(@(x) cos(4*pi*x), dom, [], [], pref);
pass(2) = abs(innerProduct(f, g)) < 10*max(f.onefun.epslevel, g.onefun.epslevel);

f = bndfun(@(x) exp(x), dom, [], [], pref);
g = bndfun(@(x) exp(-x), dom, [], [], pref);
pass(3) = abs(innerProduct(f, g) - 2) < 10*max(f.onefun.epslevel, g.onefun.epslevel);

g = bndfun(@(x) sin(x), dom, [], [], pref);
exact = exp(1)*(sin(1) - cos(1))/2 - exp(-1)*(sin(-1) - cos(-1))/2;
pass(4) = abs(innerProduct(f, g) - exact) < 10*max(f.onefun.epslevel, ...
    g.onefun.epslevel);

%%
% Check a few known properties.

f = bndfun(@(x) exp(x) - 1, dom);
g = bndfun(@(x) 1./(1 + 1i*x.^2), dom);
h = bndfun(@(x) sinh(x*exp(pi*1i/6)), dom);

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
% Check operation for array-valued bndfun objects.

f = bndfun(@(x) [sin(x) cos(x)], dom);
g = bndfun(@(x) [exp(x) 1./(1 + x.^2) airy(x)], dom);
ip = innerProduct(f, g);
exact = [0.663493666631241 0                 -0.135033172317858;
    1.933421496200713 1.365866063614065  0.592109441404267];
pass(10) = norm(ip(:) - exact(:), inf) < 10*max(f.onefun.epslevel, ...
    g.onefun.epslevel);

%%
% Check error conditions.

% Can't take the inner product of a bndfun and a non-bndfun.
try
    ip = innerProduct(f, 2); %#ok<NASGU>
    pass(11) = false;
catch ME
    pass(11) = strcmp(ME.identifier, ...
        'CHEBFUN:CHEBTECH:innerProduct:input');
end

end
