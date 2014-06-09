% Test file for singfun/innerProduct.m

function pass = test_innerProduct(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

% The order of the exponents:
a = 0.64;
b = -0.64;
c = 1.28;
d = -1.28;
p = -0.2;
q = -0.3;

% fractional pole at the left endpoint
data.exponents = [p 0];
data.singType = {'sing', 'none'};
f = singfun(@(x) (1+x).^p, data, pref);
data.exponents = [q 0];
data.singType = {'sing', 'none'};
g = singfun(@(x) (1+x).^q, data, pref);
I = innerProduct(f,g);
I_exact = 2*sqrt(2);

pass(1) = ( abs(I-I_exact) < max(get(f, 'epslevel'), get(g, 'epslevel'))*...
    abs(I_exact) );

% fractional pole at the left endpoint
data.exponents = [d 0];
data.singType = {'sing', 'none'};
f = singfun(@(x) (1+x).^d.*sin(x), data, pref);
data.exponents = [2*d 0];
data.singType = {'sing', 'none'};
g = singfun(@(x) (1+x).^(2*d), data, pref);
I = innerProduct(f,g);
pass(2) = ( isinf(I) && sign(I) == sign(sin(-1)) );

% fractional root at the right endpoint
data.exponents = [0 c];
data.singType = {'none', 'root'};
f = singfun(@(x) (1-x).^c.*cos(x), data, pref);
data.exponents = [0 a];
data.singType = {'none', 'root'};
g = singfun(@(x) (1-x).^a.*cos(x), data, pref);
I = innerProduct(f,g);
I_exact = 1.76743783779682186471;
pass(3) = ( abs(I-I_exact) < max(get(f, 'epslevel'), get(g, 'epslevel'))*...
    abs(I_exact) );

% fractional pole at the left endpoint
data.exponents = [0 b];
data.singType = {'none', 'sing'};
f = singfun(@(x) (1-x).^b.*(x.^5), data, pref);
data.exponents = [0 0];
data.singType = {'none', 'none'};
g = singfun(@(x) exp(x).*sin(5*x), data, pref);
I = innerProduct(f,g);
I_exact = -3.2185857544263774863;
pass(4) = ( abs(I-I_exact) < max(get(f, 'epslevel'), get(g, 'epslevel'))*...
    abs(I_exact) );

% a combination of fractional pole and fractional root
data.exponents = [b 0];
data.singType = {'sing', 'none'};
f = singfun(@(x) (1+x).^b.*sin(x), data, pref);
data.exponents = [0 c];
data.singType = {'none', 'root'};
g = singfun(@(x) sin(2*x).*(1-x).^c, data, pref);
I = innerProduct(f,g);
I_exact = 3.703689983503164674;
pass(5) = ( abs(I-I_exact) < max(get(f, 'epslevel'), get(g, 'epslevel'))*...
    abs(I_exact) );

% poles at different endpoints
data.exponents = [b b];
data.singType = {'sing', 'sing'};
f = singfun(@(x) sin(x).*(1-x.^2).^b, data, pref);
data.exponents = [p 0];
data.singType = {'sing', 'none'};
g = singfun(@(x) cos(x).^3.*(1+x).^p, data, pref);
I = innerProduct(f,g);
I_exact = -0.378959054771939734525;
pass(6) = ( abs(I-I_exact) < 1e1*max(get(f, 'epslevel'), get(g, 'epslevel')) );

% Check the trivial case with both vanishing alpha and beta.
data.exponents = [0 0];
data.singType = {'none', 'none'};
f = singfun(@(x) exp(x).*x.^3.*sin(2*x), data, pref);
data.exponents = [0 0];
data.singType = {'none', 'none'};
g = singfun(@(x) exp(1-x).^(3/2), data, pref);
I = innerProduct(f,g);
I_exact = 2.30589565644897950113;
pass(7) = ( abs(I-I_exact) < 1e1*max(get(f, 'epslevel'), get(g, 'epslevel'))*...
    abs(I_exact) );

% Check the complex-valued case:
f_op = @(x) (sin(x)+1i*cos(x))./((1+x).^0.4.*(1-x).^0.3);
data.exponents = [-0.4 -0.3];
data.singType = {'sing', 'sing'};
f = singfun(f_op, data, pref);
g_op = @(x) (sin(x)-1i*cos(x))./((1+x).^0.2);
data.exponents = [-0.2 0];
data.singType = {'sing', 'none'};
g = singfun(g_op, data, pref);
I = innerProduct(f,g);
I_exact = -0.66255618280005499086+0.95157967059305931745i;
pass(8) = ( abs(I-I_exact) < 1e1*max(get(f, 'epslevel'), get(g, 'epslevel'))*...
    abs(I_exact) );

end
