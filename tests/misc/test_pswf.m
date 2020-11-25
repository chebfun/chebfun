% Test file for pswf and pswfpts.m.

function pass = test_pswf(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

%% Test PSWF:
P = pswf(1:10, 4);

% Check orthogonality:
pass(1) = norm(P'*P - eye(10)) < 1e-10;

% Check number of roots:
pass(2) = numel(roots(P(:,end))) == 9;

% Check some precomputed values:
pass(3) = abs(sum(P(.5,:)) - (-0.070494617679645)) < 1e-10;

%% Test PSWFPTS:

% Check nodes are roots of P10:
[x, w] = pswfpts(9, 4);
pass(4) = norm(x - sort(roots(P(:,end)))) < 1e-10;

% CHeck integral of P(1:10):
pass(5) = norm(sum(P) - w*P(x,:)) < 1e-10;

% Check integral of cosine(k*x):
k = 51;
f = @(x) cos(k*x);
I = 2*sin(k)/k;
[x, w] = pswfpts(k, k);
err = abs(w*f(x) - I);
pass(6) = err < 1e-10;

end