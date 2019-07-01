% Test file for @chebfun/conv.m.

function pass = test_conv(pref)

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end

% Construct a few CHEBFUN objects for the tests.
f = chebfun(@(x) x, [-1 1]);
g = chebfun(@(x) sin(5*x), [2 4]);
h = chebfun(@(x) cos(2*x), [-3 1]);
p = chebfun(@(x) sin(15*x), [-1 1]);
q = chebfun(@(x) exp(cos(3*x)), [-1 1]);

%% 1. Test the commutativity.
H1 = conv(f, g);
H2 = conv(g, f);
pass(1) = norm(H1 - H2) < 100*eps;

%% 2. Test the associativity.
H1 = conv(conv(f, g), h);
H2 = conv(f, conv(g, h));
pass(2) = norm(H1 - H2) < 100*eps;

%% 3. Test the distributivity.
H1 = conv(f, (p + q));
H2 = conv(f, p) + conv(f, q);
pass(3) = norm(H1 - H2) < 100*eps;

%% 4. Test B-splines.
% Generates piecewise polynomial cardinal B-splines by a process of
% convolution, with CONV, then estimates the error compared to the exact
% solution B, using NORM. For more on generating B-splines using convolution,
% see e.g. http://en.wikipedia.org/wiki/B-spline
B = (1/6)*chebfun( {@(x) (2+x).^3, ...
    @(x)1 + 3*(1+x) + 3*(1+x).^2 - 3*(1+x).^3, ...
    @(x)1 + 3*(1-x) + 3*(1-x).^2 - 3*(1-x).^3, ...
    @(x)(2-x).^3}, -2:2 );
s = chebfun(1, [-.5 .5]);
f = s;
for k = 1:3
    f = conv(f, s);
end
pass(4) = norm(f - B) < 100*eps*vscale(f);

%% 5. Test more splines.
f = chebfun(1/2); 
g = f;
for j = 2:4
    g = conv(f, g);
end
pass(5) = abs(feval(g, 1) - 23/96) < 100*eps;

%% 6. Test the example from the HELP text:
f = chebfun(1/2); g = f;
for j = 2:4
    g = conv(f, g); 
end
g1 = g(.1);
for j = 1:4
    g = diff(g); 
end
g2 = g(.1);
err = abs(g1 - 0.332114583333333) + abs(g2 - 0);
pass(6) = err < 1e-14;
    
%% 7. Test the example from ATAP
a = 0.2885554757; b = 0.3549060246;
g = chebfun(@(x) sin(1./x).*sin(1./sin(1./x)), [a,b], 80, 'chebkind', 2);
t = 1e-7;
f = chebfun(@(x) exp(-x.^2/(4*t))/sqrt(4*pi*t),.003*[-1 1]);
h = conv(f, g);
err = abs(norm(h,2) - 0.029781437647379);
tol = 10*max(get(f, 'vscale')*eps, ...
    eps*get(g, 'vscale'));
pass(7) = err < tol;

%% An Example due to Mohsin Javed:
f = chebfun(@(x) x, [-1.5 0] );
g = chebfun(@(x) sin(x), [-1 1]);
h1 = conv(f, g);
h2 = conv(f, g, 'old');
tol = 10*max(get(f, 'vscale')*eps, ...
    eps*get(g, 'vscale'));
pass(8) = norm(h1 - h2, inf) < tol;

%% Testing Delta function convolution
% Delta funciton reproduces the function under convolution:
f = chebfun({@(x) sin(x), @(x) cos(x), @(x) sin(4*x.^2)}, [-2, -1, 0, 1] );
x = chebfun('x', [-2, 1] );
d = dirac(x);
g = conv(d, f);
g = restrict(g, [-2, 1]);
g.pointValues(1) = 2*g.pointValues(1);
g.pointValues(end) = 2*g.pointValues(end);
pass(9) = norm(f - g, inf) < tol;

% Derivative of delta function differentiates the function:
x = chebfun('x');
f = sin(x);
g = conv(f, diff(dirac(x)));
g = restrict(g, [-1, 1] );
g.pointValues(1) = 2*g.pointValues(1);
g.pointValues(end) = 2*g.pointValues(end);
pass(10) = norm(g - cos(x), inf ) < 10*tol;

% Second order ODE via delta functions and convolutions:
% g = f'' + f
f = sin(x);
g = conv(f, diff(dirac(x), 2) + dirac(x));
g = restrict(g, [-1, 1] );
g.pointValues(1) = 2*g.pointValues(1);
g.pointValues(end) = 2*g.pointValues(end);
pass(11) = norm(g , inf ) < 1e3*tol;

%% Maurice's Cox examples: 
fX1 = chebfun(@(x) exp(x), [0, log(2)]);
fX2 = chebfun(@(x) exp(x), [log(2), log(3)]);
fX3 = chebfun(@(x) exp(x), [log(3), log(4)]);
g1 = conv(fX1, fX3);
g2 = conv(g1, fX2);

g3 = conv(fX1, fX2);
g4 = conv(g3, fX3);
pass(12) = normest( g2 - g4 ) < 1e1*tol; 

g5 = conv(fX2, fX3);
g6 = conv(g5, fX1);
pass(13) = normest( g2 - g6 ) < 1e1*tol; 

%% test 'same' option
f = chebfun(@(x) exp(-x.^2), [-10 10]);
g = chebfun(@(x) exp(-x.^2), [-20 20]);
h = conv(f, f, 'same');
pass(14) = norm(h.domain([1, end]) - [-10 10], inf) < eps*10;

h = conv(f, g, 'same');
pass(15) = norm(h.domain([1, end]) - [-10 10], inf) < eps*10;

%% test quasimatrix option
f = chebfun(@sin);
g = chebfun(@cos);
fg = [f g];
gg = [g g];
h1 = conv(fg, gg);
h2 = [conv(f,g), conv(g,g)];
pass(16) = norm(h1 - h2) < eps*10;

ffg = conv([f f],g);
fgg = conv(f,[g g]);
pass(17) = norm(ffg-fgg) < 10*eps;

cheb.x;
phi = @(t) chebfun(@(x) exp(-x^2/(4*t))/sqrt(4*pi*t));
f = 1 - 0.2*x - abs(x-0.2);
fsmooth = conv(f,phi(1e-4),'same');
pass(18) = norm(f-fsmooth) < .01;

%% test small intervals
f = chebfun(@sin);
g = chebfun(@cos, [-1 1]/1e3);
h1 = conv(f, g);
h2 = conv(f,g, 'Old');
pass(19) = norm(h1 - h2) < eps*10;

f = chebfun(@sin);
g = chebfun(@cos, [-1 1]/1e6);
h1 = conv(f, g);
h2 = conv(f,g, 'Old');
pass(20) = norm(h1 - h2) < eps*10;

end

