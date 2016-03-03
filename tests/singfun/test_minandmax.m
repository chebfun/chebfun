% Test file for singfun/minandmax.m

function pass = test_minandmax(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

% The order of the exponents:
a = 0.64;
b = -0.64;
c = 1.28;
d = -1.28;

% fractional root at the left endpoint and the function value is bounded.
data.exponents = [a 0];
data.singType = {'root', 'none'};
f = singfun(@(x) (1+x).^a.*exp(x), data, pref);
[y, x] = minandmax(f);
y_exact = [0; 2^a*exp(1)];
x_exact = [-1; 1];
err_x = norm(x-x_exact, inf);
err_y = norm(y-y_exact, inf);
pass(1) = (max([err_x err_y]) < 1e1*eps*vscale(f.smoothPart));

% fractional pole at the left endpoint and the function value is unbounded.
data.exponents = [d+1 0];
data.singType = {'sing', 'none'};
f = singfun(@(x) (1+x).^d.*sin(50*pi*x), data, pref);
[y, x] = minandmax(f);
% The following exact answers are obtained using Mathematica.
y_exact = [-92.48807414703726; Inf];
x_exact = [-0.9717902218389038; -1];
err_y = y(1) - y_exact(1);
err_x = x - x_exact; 
pass(2) = (norm(err_x, inf) < 1e1*eps*vscale(f.smoothPart) &&...
    abs(err_y) < 1e4*eps*vscale(f.smoothPart)) && (y(2) == Inf);
    
% fractional root at the right endpoint and the smooth part has no roots in 
% [-1 1].
data.exponents = [0 c];
data.singType = {'none', 'root'};
f = singfun(@(x) (1-x).^c.*cos(x), data, pref);
[y, x] = minandmax(f);

% The exact maximum value and its location is found using Mathematica.
y_exact = [0; 1.511345730595468];
x_exact = [1; -0.6575681557708653];
y_err = y - y_exact;
x_err = x - x_exact;
pass(3) = (norm(y_err, inf) < 1e2*eps*vscale(f.smoothPart) &&...
    norm(x_err, inf) < 1e2*eps*vscale(f.smoothPart));

p = pref;
p.chebfuneps = max(pref.chebfuneps, 1e-14);
% no fractional pole but a root at the left endpoint.
data.exponents = [0 1+b];
data.singType = {'none', 'root'};
f = singfun(@(x) (1-x).^b.*(exp(x)-exp(1)), data, pref);
[y, x] = minandmax(f);

% The following exact answers are obtained using Mathematica.
y_exact = [-1.727141310139675; 0];
x_exact = [0.1651705232378299; 1];
y_err = y - y_exact;
x_err = x - x_exact;
pass(4) = (norm(y_err, inf) < 1e2*eps*vscale(f.smoothPart) &&...
    norm(x_err, inf) < 1e3*eps*vscale(f.smoothPart));
    
% a combination of fractional pole and fractional root.
data.exponents = [b c];
data.singType = {'sing', 'root'};
f = singfun(@(x) (1+x).^b.*sin(x).*(1-x).^c, data, pref);
[y, x] = minandmax(f);

% The exact maximum value and its location is found using Mathematica.
y_exact = [-Inf; 0.1636938399751735];
x_exact = [-1; 0.3776091222310658];
err_y = y(2) - y_exact(2);
err_x = x - x_exact; 
pass(5) = (norm(err_x, inf) < 10*eps*vscale(f.smoothPart) &&...
    abs(err_y) < 10*eps*vscale(f.smoothPart)) && (y(1) == -Inf);

end
