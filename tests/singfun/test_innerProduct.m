% Test file for singfun/innerProduct.m

function pass = test_innerProduct(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebpref();
end

% The order of the exponents:
a = 0.64;
b = -0.64;
c = 1.28;
d = -1.28;
p = -0.2;
q = -0.3;

% fractional pole at the left endpoint
f = singfun(@(x) (1+x).^p, [p 0], {'sing', 'none'}, [], [], pref);
g = singfun(@(x) (1+x).^q, [q 0], {'sing', 'none'}, [], [], pref);
I = innerProduct(f,g);
I_exact = 2*sqrt(2);

pass(1) = ( abs(I-I_exact) < max(get(f, 'epslevel'), get(g, 'epslevel'))*...
    abs(I_exact) );

% fractional pole at the left endpoint
f = singfun(@(x) (1+x).^d.*sin(x), [d 0], {'sing', 'none'}, [], [], pref);
g = singfun(@(x) (1+x).^(2*d), [2*d 0], {'sing', 'none'}, [], [], pref);
I = innerProduct(f,g);
pass(2) = ( isinf(I) && sign(I) == sign(sin(-1)) );

% fractional root at the right endpoint
f = singfun(@(x) (1-x).^c.*cos(x), [0 c], {'none', 'root'}, [], [], pref);
g = singfun(@(x) (1-x).^a.*cos(x), [0 a], {'none', 'root'}, [], [], pref);
I = innerProduct(f,g);
I_exact = 1.76743783779682186471;
pass(3) = ( abs(I-I_exact) < max(get(f, 'epslevel'), get(g, 'epslevel'))*...
    abs(I_exact) );

% fractional pole at the left endpoint
f = singfun(@(x) (1-x).^b.*(x.^5), [0 b], {'none', 'sing'}, [], [], pref);
g = singfun(@(x) exp(x).*sin(5*x), [0 0], {'none', 'none'}, [], [], pref);
I = innerProduct(f,g);
I_exact = -3.2185857544263774863;
pass(4) = ( abs(I-I_exact) < max(get(f, 'epslevel'), get(g, 'epslevel'))*...
    abs(I_exact) );

% a combination of fractional pole and fractional root
f = singfun(@(x) (1+x).^b.*sin(x), [b 0], {'sing', 'none'}, [], [], pref);
g = singfun(@(x) sin(2*x).*(1-x).^c, [0 c], {'none', 'root'}, [], [], pref);
I = innerProduct(f,g);
I_exact = 3.703689983503164674;
pass(5) = ( abs(I-I_exact) < max(get(f, 'epslevel'), get(g, 'epslevel'))*...
    abs(I_exact) );

% poles at different endpoints
f = singfun(@(x) sin(x).*(1-x.^2).^b, [b b], {'sing', 'sing'}, [], [], pref);
g = singfun(@(x) cos(x).^3.*(1+x).^p, [p 0], {'sing', 'none'}, [], [], pref);
I = innerProduct(f,g);
I_exact = -0.378959054771939734525;
pass(6) = ( abs(I-I_exact) < 1e1*max(get(f, 'epslevel'), get(g, 'epslevel')) );

% Check the trivial case with both vanishing alpha and beta.
f = singfun(@(x) exp(x).*x.^3.*sin(2*x), [0 0], {'none', 'none'}, [], [], pref);
g = singfun(@(x) exp(1-x).^(3/2), [0 0], {'none', 'none'}, [], [], pref);
I = innerProduct(f,g);
I_exact = 2.30589565644897950113;
pass(7) = ( abs(I-I_exact) < 1e1*max(get(f, 'epslevel'), get(g, 'epslevel'))*...
    abs(I_exact) );

% Check the complex-valued case:
f_op = @(x) (sin(x)+1i*cos(x))./((1+x).^0.4.*(1-x).^0.3);
f = singfun(f_op, [-0.4 -0.3], {'sing', 'sing'}, [], [], pref);
g_op = @(x) (sin(x)-1i*cos(x))./((1+x).^0.2);
g = singfun(g_op, [-0.2 0], {'sing', 'none'}, [], [], pref);
I = innerProduct(f,g);
I_exact = -0.66255618280005499086+0.95157967059305931745i;
pass(8) = ( abs(I-I_exact) < max(get(f, 'epslevel'), get(g, 'epslevel'))*...
    abs(I_exact) );

end
