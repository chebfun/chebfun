% Test file for trigtech/innerProduct.m

function pass = test_innerProduct(pref)

% Get preferences.
if ( nargin < 1 )
    pref = trigtech.techPref();
end

% Fixed arbitrary numbers to use as multiplicative constants.
alpha = -0.194758928283640 + 0.075474485412665i;
beta = -0.526634844879922 - 0.685484380523668i;

testclass = trigtech();

%%
% Spot-check a few known results.

% Orthogonality test
f = testclass.make(@(x) sin(2*pi*x), [], pref);
g = testclass.make(@(x) cos(2*pi*x), [], pref);
tol_f = 10*eps*vscale(f);
tol_g = 10*eps*vscale(g);
pass(1) = abs(innerProduct(f, g)) < max(tol_f, tol_g);

% Orthogonality test
g = testclass.make(@(x) cos(4*pi*x), [], pref);
tol_g = 10*eps*vscale(g);
pass(2) = abs(innerProduct(f, g)) < max(tol_f, tol_g);

% Easy known result
f = testclass.make(@(x) exp(cos(pi*x)), [], pref);
g = testclass.make(@(x) exp(-cos(pi*x)), [], pref);
tol_f = 10*eps*vscale(f);
tol_g = 10*eps*vscale(g);
pass(3) = abs(innerProduct(f, g) - 2) < max(tol_f, tol_g);

% Harder known result
f = testclass.make(@(x) 1 + 0*x, [], pref);
g = testclass.make(@(x) sin(pi*x).^4, [], pref);
tol_g = 10*eps*vscale(g);
exact = 3/4;
pass(4) = abs(innerProduct(f, g) - exact) < max(tol_f, tol_g);

%%
% Check a few known properties.

f = testclass.make(@(x) sin(pi*x).^4);
g = testclass.make(@(x) exp(cos(2*pi*x)));
h = testclass.make(@(x) exp(1i*4*pi*x));
tol_f = 10*eps*vscale(f);
tol_g = 10*eps*vscale(g);
tol_h = 10*eps*vscale(h);

ip1 = innerProduct(alpha*f, beta*g);
ip2 = conj(alpha)*beta*innerProduct(f, g);
pass(5) = abs(ip1 - ip2) < max(tol_f, tol_g);

ip1 = innerProduct(g, h);
ip2 = innerProduct(h, g);
pass(6) = abs(ip1 - conj(ip2)) < max(tol_g, tol_h);

ip1 = innerProduct(f + g, h);
ip2 = innerProduct(f, h) + innerProduct(g, h);
pass(7) = abs(ip1 - ip2) < max([tol_f tol_g tol_h]);

ip1 = innerProduct(f, g + h);
ip2 = innerProduct(f, g) + innerProduct(f, h);
pass(8) = abs(ip1 - ip2) < max([tol_f, tol_g, tol_h]);

nf2 = innerProduct(f, f);
ng2 = innerProduct(g, g);
nh2 = innerProduct(h, h);
n2vals = [nf2 ; ng2 ; nh2];
pass(9) = isreal(n2vals) && all(n2vals >= 0);

%% 
% Check operation for array-valued trigtech objects.

f = testclass.make(@(x) [cos(pi*sin(2*pi*x)) exp(sin(2*pi*x)) sin(pi*sin(pi*x))]);
g = testclass.make(@(x) [exp(cos(pi*x)) exp(-sin(2*pi*x)) cos(pi*x)]);
tol_f = 10*max(vscale(f)*eps);
tol_g = 10*max(vscale(g)*eps);
ip = innerProduct(f, g);
exact = [-0.765066454912991 -1.032310819204998 0;
          3.204359383938947  2                 0;
          0                  0                 0];
pass(10) = norm(ip(:) - exact(:), inf) < max(tol_f, tol_g);

%%
% Check error conditions.

% Can't take the inner product of a trigtech and a non-chebtech.
try
    ip = innerProduct(f, 2); %#ok<NASGU>
    pass(11) = false;
catch ME
    pass(11) = strcmp(ME.identifier, ...
        'CHEBFUN:TRIGTECH:innerProduct:input');
end

end
