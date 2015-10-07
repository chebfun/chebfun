function pass = test_inv(pref)
% This test constructs two CHEBFUN objects and uses INV() to invert them. It
% checks the inverse calculated is accurate.

% Taken from Chebfun v4 test, invtest.m, by Nick Hale  07/06/2009.

if ( nargin == 0 ) 
    pref = chebfunpref();
end

algoList = {'roots', 'newton', 'bisection', 'regulafalsi', 'illinois', 'brent'};

for k = 1:6
    x = chebfun('x');
    f = sin(x);
    g = chebfun(@(x) asin(x), [sin(-1), sin(1)]);
    f_inv = inv(f, pref, 'algorithm', algoList{k});
    tol = 100*eps*vscale(f_inv);
    pass(k,1) = norm(g - f_inv, inf) < tol;

    g = inv(f_inv, 'algorithm', algoList{k});
    [f, g] = tweakDomain(f, g, 2e-15);
    pass(k,2) = norm(f - g, inf) < tol;

    x = chebfun('x');
    f = chebfun(@(x) sausagemap(x));
    f_inv = inv(f, pref, 'algorithm', algoList{k});
    tol = 100*eps*vscale(f_inv);
    pass(k,3) = norm(f(f_inv) - x, inf) + norm(f_inv(f) - x, inf) < tol;

    % Check the 'monocheck' and 'rangecheck' options.
    f = chebfun(@exp);
    f_inv = inv(f, pref, 'algorithm', algoList{k}, 'monocheck', 'on', ...
        'rangecheck', 'on');
    tol = 100*eps*vscale(f_inv);
    xx = linspace(f_inv.domain(1), f_inv.domain(end), 20);
    pass(k,4) = norm(f(f_inv(xx)) - xx, inf) < tol;
    pass(k,5) = all(abs(f_inv.domain([1, end]) - exp([-1 1])) < 10*eps);
end    

% Check that 'monocheck' fails for a non-monotonic function.
try
    f = chebfun(@(x) x.^2);
    f_inv = inv(f, pref, 'splitting', 'on', 'monocheck', 'on');
    pass(:,6) = false;
catch ME
    if ( strcmpi(ME.identifier, ...
            'CHEBFUN:CHEBFUN:inv:doMonoCheck:notMonotonic') )
        pass(:,6) = true;
    else
        pass(:,6) = false;
    end
end

% Test the example from the help text:
x = chebfun('x');
f = x + .5*abs(x) + .6*sign(x-.5);
g = inv(f);
xx = linspace(f.domain(1), f.domain(end), 10);
pass(:,7) = norm(g(f(xx)) - xx, inf) < 10*vscale(g)*eps;

% Test the inverse of a decreasing function (see #1098).
x = chebfun('x');
f = -sin(x);
g = inv(f);
xx = linspace(f.domain(1), f.domain(end), 10);
pass(:,8) = norm(g(f(xx)) - xx, inf) < 10*vscale(g)*eps;

end

function g = sausagemap(s,d)
if ( nargin < 2 )
    d = 9; % This can be adjusted
end 
c = zeros(1,d+1);
c(d:-2:1) = [1 cumprod(1:2:d-2)./cumprod(2:2:d-1)]./(1:2:d);
c = c/sum(c); g = polyval(c,s);
cp = c(1:d).*(d:-1:1);
end
