% Test file for singfun constructor.

function pass = test_singfun_constructor(pref)

% Get preferences:
if ( nargin < 1 )
    pref = chebfunpref();
end

%%
% Select some random points as sample points
% These random points are in [-1, 1]
seedRNG(7890)
x = -1 + 2*rand(1, 100);
x = sort(x);

% Some arbitrary values to use for exponents.
a = 0.338745372057174;
b = 0.561224728136042;
a_int = 2;
b_int = 5;


%% Test calling syntax when the user provides exponents

% Negative fractional exponents
fh = @(x) sin(x)./((1+x).^a.*(1-x).^b);
data = struct();
data.exponents = [-a, -b];
f = singfun(fh, data, pref);
data.exponents = [-a, -b];
data.singType = {'sing', 'sing'};
g = singfun(fh, data, pref);
pass(1) = isequal(f,g);
pass(2) = ~any(f.exponents + [a,b]);
pass(3) = ~any(g.exponents + [a,b]);
pass(4) = norm(feval(fh,x) - feval(f,x), inf) < 1e2*eps;

%%
% Positive fractional exponents
fh = @(x) sin(x).*(1+x).^a.*(1-x).^b;
data = struct();
data.exponents = [a, b];
f = singfun(fh, data, pref);
data.exponents = [a, b];
data.singType = {'root', 'root'};
g = singfun(fh, data, pref);
pass(5) = isequal(f,g);
pass(6) = ~any(f.exponents - [a,b]);
pass(7) = ~any(g.exponents - [a,b]);
pass(8) = norm(feval(fh,x) - feval(f,x), inf) < 1e1*eps;

%%
% Negative integer exponents
fh = @(x) exp(x)./((1+x).^a_int.*(1-x).^b_int);
data = struct();
data.exponents = [-a_int, -b_int];
f = singfun(fh, data, pref);
data.exponents = [-a_int, -b_int];
data.singType = {'pole', 'pole'};
g = singfun(fh, data, pref);
pass(9) = isequal(f,g);
pass(10) = ~any(f.exponents + [a_int, b_int]);
pass(11) = ~any(g.exponents + [a_int, b_int]);
% don't check near end-points
xx = x(20:80);
pass(12) = norm(feval(fh,xx) - feval(f,xx), inf) < 1e2*eps;

%% Test Syntax and construction when the user doesn't provide exponents
%
% Negative fractional exponents
fh = @(x) exp(sin(x))./((1+x).^a.*(1-x).^b);
f = singfun(fh);
pass(13) = norm(f.exponents + [a,b], inf) < pref.blowupPrefs.exponentTol;
pass(14) = norm(feval(fh,x) - feval(f,x), inf) < 1e5*eps;
    
%%
% Positive fractional exponents
fh = @(x) sin(exp(cos(x))).*(1+x).^a.*(1-x).^b;
f = singfun(fh);
pass(15) = norm(f.exponents - [a,b], inf) < pref.blowupPrefs.exponentTol;
pass(16) = norm(feval(fh,x) - feval(f,x), inf) < 1e4*eps;
    
%%
% Negative integer exponents
fh = @(x) exp(sin(x.^2))./((1+x).^a_int.*(1-x).^b_int);
f = singfun(fh);
pass(17) = norm(f.exponents + [a_int, b_int], inf) < pref.blowupPrefs.exponentTol;
xx = x(20:80);
pass(18) = norm(feval(fh,xx) - feval(f,xx), inf) < 1e3*eps;


%%
% Construction with smoothfuns:
f = smoothfun.constructor( @(x) sin(x));
s = singfun(f);
pass(19) = iszero(f - s.smoothPart);
data = struct();
data.exponents = [-1.5, -1];
data.singType = {'sing', 'sing'};
s = singfun(f, data, pref);
pass(20) = iszero(f - s.smoothPart);
pass(21) = norm(s.exponents - [-1.5, -1], inf) < pref.blowupPrefs.exponentTol;

%%
% Construction from double:
f = singfun(42);
pass(22) = iszero(f - 42);
data = struct();
data.exponents = [1.5, 1];
data.singType = {'sing', 'sing'};
f = singfun(42, data, pref);
g = singfun(@(x) 42+0*x);
g.exponents = [1.5, 1];
pass(23) = iszero(f - g);
pass(24) = norm(s.exponents - [-1.5, -1], inf) < pref.blowupPrefs.exponentTol;

end
