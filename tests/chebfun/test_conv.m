function pass = test_conv(pref)

% Grab some preferences
if ( nargin == 0 )
    pref = chebpref();
end

% Construct a few CHEBFUN objects for the tests.
f = chebfun('x', [-1 1]);
g = chebfun('sin(5*x)', [2 4]);
h = chebfun('cos(2*x)', [-3 1]);
p = chebfun('sin(15*x)', [-1 1]);
q = chebfun('exp(cos(3*x))', [-1 1]);

%% 1. test the commutativity
H1 = conv(f, g);
H2 = conv(g, f);
pass(1) = norm(H1 - H2) < 100*epslevel(H1);

%% 2. test the associativity
H1 = conv(conv(f, g), h);
H2 = conv(f, conv(g, h));
pass(2) = norm(H1 - H2) < 100*epslevel(H1);

%% 3. test the distributivity
H1 = conv(f, (p + q));
H2 = conv(f, p) + conv(f, q);
pass(3) = norm(H1 - H2) < 100*epslevel(H1);

%% 4. test B-splines 
% Generates piecewise polynomial cardinal B-splines by a process of convolution,
% with CONV, then estimates the error compared to the exact solution B, using
% Chebfun's NORM command. For more on generating B-splines using convolution,
% see e.g. http://en.wikipedia.org/wiki/B-spline
B = (1/6)*chebfun( {@(x) (2+x).^3, ...
    @(x)1 + 3*(1+x) + 3*(1+x).^2 - 3*(1+x).^3, ...
    @(x)1 + 3*(1-x) + 3*(1-x).^2 - 3*(1-x).^3, ...
    @(x)(2-x).^3}, -2:2 );
s = chebfun(1, [-.5 .5]);
f = s;
for k = 1:3, f = conv(f, s); end
pass(4) = norm( f - B ) < 100*epslevel(f)*vscale(f);

%% 5. test more splines
f = chebfun(1/2); 
g = f;
for j = 2:4, g = conv(f, g); end
pass(5) = numel(g.domain) == 5 && all(g.domain == -4:2:4) && ...
    abs(feval(g, 1) - 23/96) < 100*epslevel(g);

end


