% Test file for singfun/sum.m

function pass = test_sum(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

% The order of the exponents:
a = 0.64;
b = -0.64;
c = 1.28;
d = -1.28;

% fractional root at the left endpoint
data.exponents = [a 0];
data.singType = {'root', 'none'};
f = singfun(@(x) (1+x).^a.*exp(x), data, pref);
I = sum(f);
I_exact = 2.7263886326217359442;

pass(1) = ( abs(I-I_exact) < get(f, 'epslevel')*abs(I_exact) );

% fractional pole at the left endpoint
data.exponents = [d 0];
data.singType = {'sing', 'none'};
f = singfun(@(x) (1+x).^d.*sin(x), data, pref);
I = sum(f);
pass(2) = ( isinf(I) && sign(I) == sign(sin(-1)) );

% fractional root at the right endpoint
data.exponents = [0 c];
data.singType = {'none', 'root'};
f = singfun(@(x) (1-x).^c.*cos(x), data, pref);
I = sum(f);
I_exact = 1.7756234306626192717;
pass(3) = ( abs(I-I_exact) < get(f, 'epslevel')*abs(I_exact) );

% fractional pole at the left endpoint
data.exponents = [0 b];
data.singType = {'none', 'sing'};
f = singfun(@(x) (1-x).^b.*(x.^5), data, pref);
I = sum(f);
I_exact = 1.2101935306745953520;
pass(4) = ( abs(I-I_exact) < get(f, 'epslevel')*abs(I_exact) );

% a combination of fractional pole and fractional root
data.exponents = [b c];
data.singType = {'sing', 'root'};
f = singfun(@(x) (1+x).^b.*sin(x).*(1-x).^c, data, pref);
I = sum(f);
I_exact = -3.8210796477539148513;
pass(5) = ( abs(I-I_exact) < get(f, 'epslevel')*abs(I_exact) );

% test the case that pole orders at the ends are same
data.exponents = [b b];
data.singType = {'sing', 'sing'};
f = singfun(@(x) sin(x).*(1-x.^2).^b, data, pref);
I = sum(f);
I_exact = 0;
pass(6) = ( abs(I-I_exact) < get(f, 'epslevel') );

% Check the trivial case with both vanishing alpha and beta.
data.exponents = [0 0];
data.singType = {'none', 'none'};
f = singfun(@(x) exp(x).*x.^3.*sin(2*x), data, pref);
I = sum(f);
I_exact = 0.644107794617991224;
err = abs(I-I_exact);
tol = 10*get(f, 'epslevel')*abs(I_exact);
pass(7) = ( err < tol );

end
