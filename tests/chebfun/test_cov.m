% Test file for @chebfun/cov.m.

function pass = test_cov(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

% Check empty case.
pass(1) = isnan(cov(chebfun()));

% Check cov() for a single chebfun argument.
f = chebfun({@(x) exp(4*pi*1i*x), @exp, @exp}, [-1 0 0.5 1], pref);
pass(2) = isequal(cov(f), var(f));

ft = f.';
pass(3) = isequal(cov(ft), var(ft));

f = chebfun(@(x) [sin(x) cos(x) exp(x)], [-1 -0.5 0 0.5 1], pref);
css = (1 - sin(2)/2)/2;
csc = 0;
cse = (sin(1) + cos(1) + exp(2)*(sin(1) - cos(1)))/(4*exp(1));
ccc = ((1 + sin(1)*cos(1))/2 - sin(1)^2);
cce = (sin(1) - cos(1) + exp(2)*(sin(1) + cos(1)))/(4*exp(1)) ...
    - sin(1)*(exp(1) - exp(-1))/2;
cee = (1 - exp(-2))/2;
covf_exact = [css csc cse
              csc ccc cce
              cse cce cee];
covf = cov(f);
pass(4) = norm(covf(:) - covf_exact(:)) < 10*vscale(f)*eps;

ft = f.';
covft = cov(ft);
covft_exact = covf_exact.';
pass(5) = norm(covft(:) - covft_exact(:)) < 10*vscale(ft)*eps;

% Check cov() for two chebfun arguments.
covff = cov(f, f);
pass(6) = norm(covff(:) - covf_exact(:)) < 10*vscale(f)*eps;

f = chebfun(@(x) [x sin(x)], pref);
g = chebfun(@(x) [exp(x) exp(2*pi*1i*x)], [-1 -0.5 0 0.5 1], pref);
covfg_exact = [1/exp(1) -1/(2*pi*1i) ; cse -2*pi*1i*sin(1)/(1 - 4*pi^2)];
covfg = cov(f, g);
pass(7) = norm(covfg(:) - covfg_exact(:)) < ...
    10*max(vscale(f)*eps, vscale(g)*eps);

ft = f.';
gt = g.';
covfgt_exact = covfg_exact.';
covfgt = cov(ft, gt);
pass(8) = norm(covfgt(:) - covfgt_exact(:)) < ...
    10*max(vscale(f)*eps, vscale(g)*eps);

% Check error conditions.
try
    f = chebfun(@sin, pref);
    g = chebfun(@(x) [sin(x) cos(x)], pref);
    C = cov(f, g);
    pass(9) = false;
catch ME
    pass(9) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:cov:size');
end

try
    C = cov(f, g, 1);
    pass(10) = false;
catch ME
    pass(10) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:cov:nargin');
end

end
