% Test file for singfun/roots.m

function pass = test_roots(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

% The order of the exponents:
a = 0.56;
b = -0.56;
c = 1.28;
d = -1.28;

% fractional root at the left endpoint and the smooth part has no roots in 
% [-1 1].
data.exponents = [a 0];
data.singType = {'root', 'none'};
f = singfun(@(x) (1+x).^a.*exp(x), data, pref);
r = roots(f);
r_exact = -1;
err = r - r_exact;
pass(1) = (norm(err, inf) < 10*get(f, 'vscale')*get(f, 'epslevel'));

% fractional pole at the left endpoint and note that the left endpoint is 
% not a root, even though the smooth part has a root there.
data.exponents = [d+1 0];
data.singType = {'sing', 'none'};
f = singfun(@(x) (1+x).^d.*sin(50*pi*x), data, pref);
r = roots(f);
r_exact = (-1+1/50:(1/50):1).';
err = r - r_exact;
pass(2) = (norm(err, inf) < 10*get(f, 'vscale')*get(f, 'epslevel'));

% fractional root at the right endpoint and the smooth part has no roots in [-1
% 1].
data.exponents = [0 c];
data.singType = {'none', 'sing'};
f = singfun(@(x) (1-x).^c.*cos(x), data, pref);
r = roots(f);
r_exact = 1;
err = r - r_exact;
pass(3) = (norm(err, inf) < 10*get(f, 'vscale')*get(f, 'epslevel'));

% no fractional pole but a root at the right endpoint.
data.exponents = [0 1+b];
data.singType = {'none', 'root'};
f = singfun(@(x) (1-x).^b.*(exp(x)-exp(1)), data, pref);
r = roots(f);
r_exact = 1;
err = r - r_exact;
pass(4) = (norm(err, inf) < 10*get(f, 'vscale')*get(f, 'epslevel'));

% a combination of fractional pole and fractional root.
data.exponents = [b c];
data.singType = {'sing', 'root'};
f = singfun(@(x) (1+x).^b.*sin(x).*(1-x).^c, data, pref);
r = roots(f);
r_exact = [0; 1];
err = r - r_exact;
pass(5) = (norm(err, inf) < 10*get(f, 'vscale')*get(f, 'epslevel'));

% Check the case with roots close to endpoints.
p = 1-1e-14;
data.exponents = [b b];
data.singType = {'sing', 'sing'};
f = singfun(@(x) (1+x).^b.*sin(x-p).*(1-x).^b, data, pref);
r = roots(f);
r_exact = p;
err = r - r_exact;
pass(6) = (norm(err, inf) < 10*get(f, 'vscale')*get(f, 'epslevel'));

end
